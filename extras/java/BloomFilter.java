/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.1
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package btllib;

public class BloomFilter {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected BloomFilter(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(BloomFilter obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  @SuppressWarnings("deprecation")
  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        btllibJNI.delete_BloomFilter(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public BloomFilter(long size) {
    this(btllibJNI.new_BloomFilter(size), true);
  }

  public void insert(SWIGTYPE_p_std__vectorT_uint64_t_t hashes) {
    btllibJNI.BloomFilter_insert__SWIG_0(swigCPtr, this, SWIGTYPE_p_std__vectorT_uint64_t_t.getCPtr(hashes));
  }

  public void insert(SWIGTYPE_p_uint64_t hashes) {
    btllibJNI.BloomFilter_insert__SWIG_1(swigCPtr, this, SWIGTYPE_p_uint64_t.getCPtr(hashes));
  }

  public boolean contains(SWIGTYPE_p_std__vectorT_uint64_t_t hashes) {
    return btllibJNI.BloomFilter_contains__SWIG_0(swigCPtr, this, SWIGTYPE_p_std__vectorT_uint64_t_t.getCPtr(hashes));
  }

  public boolean contains(SWIGTYPE_p_uint64_t hashes) {
    return btllibJNI.BloomFilter_contains__SWIG_1(swigCPtr, this, SWIGTYPE_p_uint64_t.getCPtr(hashes));
  }

}
