#ifndef BTLLIB_NUM_QUEUE_HPP
#define BTLLIB_NUM_QUEUE_HPP

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace btllib {

const size_t NUM_QUEUE_BUFFER_SIZE = 16384;

// Surrounds pieces of data in the buffer with a busy mutex
// for exclusive access
template<typename T>
struct Slot
{
	T data;
	std::mutex busy;
	bool occupied = false;
	std::condition_variable occupancyChanged;
	long lastTenant = -1; // Required to ensure read order
};

template<typename T>
class NumQueue
{

  public:
	size_t elements() const { return element_count; }

	void close()
	{
		closed = true;
		for (Slot<T>& slot : this->slots) {
			slot.occupancyChanged.notify_all();
		}
	}

	bool is_closed() const { return closed; }

  protected:
	std::vector<Slot<T>> slots{ NUM_QUEUE_BUFFER_SIZE };
	size_t read_counter = 0;
	std::atomic<size_t> element_count{ 0 };
	std::atomic<bool> closed{ false };
};

template<typename T>
class InputNumQueue : public NumQueue<T>
{

  public:

	void write(T& data)
	{
		size_t num = data.num;
		Slot<T>& target = this->slots[num % NUM_QUEUE_BUFFER_SIZE];
		std::unique_lock<std::mutex> busy_lock(target.busy);
		target.occupancyChanged.wait(busy_lock, [&] { return !target.occupied; });
		target.data = std::move(data);
		target.occupied = true;
		target.occupancyChanged.notify_one();
		++(this->element_count);
	}

	void read(T& data)
	{
		static std::mutex read_mutex;
		std::unique_lock<std::mutex> read_lock(read_mutex);

		Slot<T>& target = this->slots[this->read_counter % NUM_QUEUE_BUFFER_SIZE];
		std::unique_lock<std::mutex> busy_lock(target.busy);
		target.occupancyChanged.wait(busy_lock, [&] { return target.occupied || this->closed; });
		if (this->closed) {
			return;
		}
		++(this->read_counter);

		read_lock.unlock();

		data = std::move(target.data);
		target.occupied = false;
		target.occupancyChanged.notify_one();
		--(this->element_count);
	}
};

template<typename T>
class OutputNumQueue : public NumQueue<T>
{

  public:
	void write(T& data)
	{
		size_t num = data.num;
		Slot<T>& target = this->slots[num % NUM_QUEUE_BUFFER_SIZE];
		std::unique_lock<std::mutex> busy_lock(target.busy);
		target.occupancyChanged.wait(
		    busy_lock, [&] { return !target.occupied && (num - target.lastTenant <= NUM_QUEUE_BUFFER_SIZE); });
		target.data = std::move(data);
		target.occupied = true;
		target.lastTenant = num;
		target.occupancyChanged.notify_all();
		++(this->element_count);
	}

	void read(T& data)
	{
		Slot<T>& target = this->slots[this->read_counter % NUM_QUEUE_BUFFER_SIZE];
		std::unique_lock<std::mutex> busy_lock(target.busy);
		target.occupancyChanged.wait(busy_lock, [&] { return target.occupied || this->closed; });
		if (this->closed) {
			return;
		}
		++(this->read_counter);
		data = std::move(target.data);
		target.occupied = false;
		target.occupancyChanged.notify_all();
		--(this->element_count);
	}
};

} // namespace btllib

#endif