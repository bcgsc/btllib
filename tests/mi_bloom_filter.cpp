#include "btllib/mi_bloom_filter.hpp"
#include "btllib/nthash.hpp"

#include "helpers.hpp"

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>
#include <stdlib.h>     /* srand, rand */
int
main()
{
  std::cerr << "Testing multi-indexed BloomFilter" << std::endl;
  btllib::MIBloomFilter<uint8_t> mi_bf_1(1024 * 1024, 3, "ntHash");
  mi_bf_1.insert_bv({ 1, 10, 100 });
  mi_bf_1.insert_bv({ 100, 200, 300 });
  mi_bf_1.complete_bv_insertion();

  TEST_ASSERT(mi_bf_1.bv_contains({ 1, 10, 100 }));
  TEST_ASSERT(mi_bf_1.bv_contains({ 100, 200, 300 }));
  TEST_ASSERT(!mi_bf_1.bv_contains({ 1, 20, 100 }));

  uint8_t ID_1 = 12;
  mi_bf_1.insert_id({ 1, 10, 100 }, ID_1);
  
  std::vector<uint8_t> results_1 = mi_bf_1.get_id({ 1, 10, 100 });
  for(auto& id : results_1) {
  	TEST_ASSERT_EQ(id, ID_1);
  }

  std::cerr << "multi-indexed BloomFilter ID count test" << std::endl;
  std::cerr << "Testing ID counting" << std::endl;
  bool include_saturated = true;

  TEST_ASSERT(mi_bf_1.get_ID_occurence_count(include_saturated)[ID_1] == 3)


  std::cerr << "Testing multi-indexed BloomFilter random sampling" << std::endl;  
  int dna_length = 100000;
  std::string random_dna = get_random_seq(dna_length);
  int expected_id_count = dna_length / 4;
  double tolerance = 0.1;

  int counter = 0;
  btllib::MIBloomFilter<uint8_t> mi_bf_2(256 * 1024 * 1024, 1, "ntHash");
  for(btllib::NtHash nthash(random_dna, 1, 15); nthash.roll(); counter++){
	mi_bf_2.insert_bv(nthash.hashes());
  }
  mi_bf_2.complete_bv_insertion();
  uint8_t ID_array[4] = {0, 1, 2, 3};
  for(auto& id : ID_array){
  	uint8_t ID = id;
	for(btllib::NtHash nthash(random_dna, 1, 15); nthash.roll(); counter++){
		mi_bf_2.insert_id(nthash.hashes(), ID);
  	}
  }
  std::vector<uint8_t> results_2(1);
  std::vector<uint32_t> total_counter(4, 0);
  for(btllib::NtHash nthash(random_dna, 1, 15); nthash.roll(); counter++){
  	results_2 = mi_bf_2.get_id(nthash.hashes());
  	for(auto& res : results_2){
	      	total_counter[res]++;
	}
  }		
  for(auto& count : total_counter){
	TEST_ASSERT(count < expected_id_count + (expected_id_count * tolerance));
	TEST_ASSERT(count > expected_id_count - (expected_id_count * tolerance));
  }
  TEST_ASSERT(mi_bf_2.get_pop_saturated_cnt() == 0); // testing no saturation
  
 

  std::cerr << "Testing multi-indexed BloomFilter saving." << std::endl; 
  mi_bf_2.save("test.mibf");

  std::cerr << "Testing multi-indexed BloomFilter reading." << std::endl;
  btllib::MIBloomFilter<uint8_t> mi_bf_3("test.mibf");
  
  std::vector<uint32_t> total_counter_2(4, 0);
  for(btllib::NtHash nthash(random_dna, 1, 15); nthash.roll(); counter++){
        results_2 = mi_bf_3.get_id(nthash.hashes());
        for(auto& res : results_2){
                total_counter_2[res]++;
        }
  }
  
  for(uint k=0; k < total_counter_2.size(); k++){
	TEST_ASSERT(total_counter_2[k] == total_counter[k]);
  }
  // Test mi-Bf is still insertable.
  uint8_t ID = 3;
  for(btllib::NtHash nthash(random_dna, 1, 15); nthash.roll(); counter++){
  	mi_bf_3.insert_id(nthash.hashes(), ID);
  }

  std::cerr << "Testing multi-indexed BloomFilter saturation" << std::endl;
  btllib::MIBloomFilter<uint8_t> mi_bf_4(1024 * 1024, 3, "ntHash");
  mi_bf_4.insert_bv({ 1, 10, 100 });
  mi_bf_4.insert_bv({ 100, 200, 300 });
  mi_bf_4.insert_bv({ 500, 1000, 2000 });
  mi_bf_4.complete_bv_insertion();

  ID_1 = 1;
  uint8_t ID_2 = 2, ID_3 = 3, ID_4 = 4;
  mi_bf_4.insert_id({ 1, 10, 100 }, ID_1);
  mi_bf_4.insert_id({ 1, 10, 100 }, ID_2);
  mi_bf_4.insert_id({ 500, 1000, 2000 }, ID_1);
  mi_bf_4.insert_id({ 500, 1000, 2000 }, ID_2);
  mi_bf_4.insert_id({ 500, 1000, 2000 }, ID_3);
  mi_bf_4.insert_id({ 500, 1000, 2000 }, ID_4);
  mi_bf_4.complete_id_insertion();

  mi_bf_4.insert_saturation({ 1, 10, 100 }, ID_1);
  mi_bf_4.insert_saturation({ 1, 10, 100 }, ID_2);
  mi_bf_4.insert_saturation({ 500, 1000, 2000 }, ID_1);
  mi_bf_4.insert_saturation({ 500, 1000, 2000 }, ID_2);
  mi_bf_4.insert_saturation({ 500, 1000, 2000 }, ID_3);
  mi_bf_4.insert_saturation({ 500, 1000, 2000 }, ID_4);
  
  // both should be represented unsaturated
  bool ID_1_found = false, ID_2_found = false, ID_3_found = false, ID_4_found = false;
  std::vector<uint8_t> results_3;
  results_3 =  mi_bf_4.get_id({ 1, 10, 100 });
  
  // both of ID's should be found.
  for(auto& id : results_3){
  	ID_1_found = id == ID_1 ? true : ID_1_found;
	ID_2_found = id == ID_1 ? true : ID_2_found;
  }

  TEST_ASSERT(ID_1_found);
  TEST_ASSERT(ID_2_found);

  // all should be saturated
  ID_1_found = false, ID_2_found = false;
  results_3 =  mi_bf_4.get_id({ 500, 1000, 2000 });

  for(auto& id : results_3){
	if(id < mi_bf_4.MASK){continue;}
        ID_1_found = (id & mi_bf_4.ANTI_MASK) == ID_1 ? true : ID_1_found;
        ID_2_found = (id & mi_bf_4.ANTI_MASK) == ID_2 ? true : ID_2_found;
	ID_3_found = (id & mi_bf_4.ANTI_MASK) == ID_3 ? true : ID_3_found;
	ID_4_found = (id & mi_bf_4.ANTI_MASK) == ID_4 ? true : ID_4_found;
  }
  //one must be absent others must be saturated.
  TEST_ASSERT(
		(!ID_1_found && (ID_2_found && ID_3_found && ID_4_found))
		|| (!ID_2_found && (ID_3_found && ID_4_found && ID_1_found))
		|| (!ID_3_found && (ID_4_found && ID_1_found && ID_2_found))
		|| (!ID_4_found && (ID_1_found && ID_2_found && ID_3_found))
  );
  // get pop saturated count should return positive integer.
  TEST_ASSERT(mi_bf_4.get_pop_saturated_cnt() > 0);

  return 0;
}
