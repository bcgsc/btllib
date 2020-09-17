/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package btllib;

public class btllib {
  public static SWIGTYPE_p_uint8_t getCP_OFF() {
    return new SWIGTYPE_p_uint8_t(btllibJNI.CP_OFF_get(), true);
  }

  public static int getMULTISHIFT() {
    return btllibJNI.MULTISHIFT_get();
  }

  public static SWIGTYPE_p_uint64_t getMULTISEED() {
    return new SWIGTYPE_p_uint64_t(btllibJNI.MULTISEED_get(), true);
  }

  public static SWIGTYPE_p_uint64_t getSEED_A() {
    return new SWIGTYPE_p_uint64_t(btllibJNI.SEED_A_get(), true);
  }

  public static SWIGTYPE_p_uint64_t getSEED_C() {
    return new SWIGTYPE_p_uint64_t(btllibJNI.SEED_C_get(), true);
  }

  public static SWIGTYPE_p_uint64_t getSEED_G() {
    return new SWIGTYPE_p_uint64_t(btllibJNI.SEED_G_get(), true);
  }

  public static SWIGTYPE_p_uint64_t getSEED_T() {
    return new SWIGTYPE_p_uint64_t(btllibJNI.SEED_T_get(), true);
  }

  public static SWIGTYPE_p_uint64_t getSEED_N() {
    return new SWIGTYPE_p_uint64_t(btllibJNI.SEED_N_get(), true);
  }

  public static int getASCII_SIZE() {
    return btllibJNI.ASCII_SIZE_get();
  }

  public static SWIGTYPE_p_uint64_t getSeed_tab() {
    long cPtr = btllibJNI.seed_tab_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t getA33R() {
    long cPtr = btllibJNI.A33R_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t getA31L() {
    long cPtr = btllibJNI.A31L_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t getC33R() {
    long cPtr = btllibJNI.C33R_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t getC31L() {
    long cPtr = btllibJNI.C31L_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t getG33R() {
    long cPtr = btllibJNI.G33R_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t getG31L() {
    long cPtr = btllibJNI.G31L_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t getT33R() {
    long cPtr = btllibJNI.T33R_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t getT31L() {
    long cPtr = btllibJNI.T31L_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t getN33R() {
    long cPtr = btllibJNI.N33R_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t getN31L() {
    long cPtr = btllibJNI.N31L_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static void setMs_tab_33r(SWIGTYPE_p_p_uint64_t value) {
    btllibJNI.ms_tab_33r_set(SWIGTYPE_p_p_uint64_t.getCPtr(value));
  }

  public static SWIGTYPE_p_p_uint64_t getMs_tab_33r() {
    long cPtr = btllibJNI.ms_tab_33r_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_p_uint64_t(cPtr, false);
  }

  public static void setMs_tab_31l(SWIGTYPE_p_p_uint64_t value) {
    btllibJNI.ms_tab_31l_set(SWIGTYPE_p_p_uint64_t.getCPtr(value));
  }

  public static SWIGTYPE_p_p_uint64_t getMs_tab_31l() {
    long cPtr = btllibJNI.ms_tab_31l_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint8_t getRC_CONVERT_TAB() {
    long cPtr = btllibJNI.RC_CONVERT_TAB_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint8_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint8_t getCONVERT_TAB() {
    long cPtr = btllibJNI.CONVERT_TAB_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint8_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t getDIMER_TAB() {
    long cPtr = btllibJNI.DIMER_TAB_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t getTRIMER_TAB() {
    long cPtr = btllibJNI.TRIMER_TAB_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t getTETRAMER_TAB() {
    long cPtr = btllibJNI.TETRAMER_TAB_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_uint64_t(cPtr, false);
  }

  public static SWIGTYPE_p_uint64_t rol1(SWIGTYPE_p_uint64_t v) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.rol1(SWIGTYPE_p_uint64_t.getCPtr(v)), true);
  }

  public static SWIGTYPE_p_uint64_t rolx(SWIGTYPE_p_uint64_t v, long x) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.rolx(SWIGTYPE_p_uint64_t.getCPtr(v), x), true);
  }

  public static SWIGTYPE_p_uint64_t ror1(SWIGTYPE_p_uint64_t v) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.ror1(SWIGTYPE_p_uint64_t.getCPtr(v)), true);
  }

  public static SWIGTYPE_p_uint64_t rol31(SWIGTYPE_p_uint64_t v, long s) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.rol31(SWIGTYPE_p_uint64_t.getCPtr(v), s), true);
  }

  public static SWIGTYPE_p_uint64_t rol33(SWIGTYPE_p_uint64_t v, long s) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.rol33(SWIGTYPE_p_uint64_t.getCPtr(v), s), true);
  }

  public static SWIGTYPE_p_uint64_t swapbits033(SWIGTYPE_p_uint64_t v) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.swapbits033(SWIGTYPE_p_uint64_t.getCPtr(v)), true);
  }

  public static SWIGTYPE_p_uint64_t swapbits3263(SWIGTYPE_p_uint64_t v) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.swapbits3263(SWIGTYPE_p_uint64_t.getCPtr(v)), true);
  }

  public static SWIGTYPE_p_uint64_t swapxbits033(SWIGTYPE_p_uint64_t v, long x) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.swapxbits033(SWIGTYPE_p_uint64_t.getCPtr(v), x), true);
  }

  public static SWIGTYPE_p_uint64_t ntf64(String kmer_seq, long k) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.ntf64__SWIG_0(kmer_seq, k), true);
  }

  public static SWIGTYPE_p_uint64_t ntr64(String kmer_seq, long k) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.ntr64__SWIG_0(kmer_seq, k), true);
  }

  public static SWIGTYPE_p_uint64_t ntf64(SWIGTYPE_p_uint64_t fh_val, long k, short char_out, short char_in) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.ntf64__SWIG_1(SWIGTYPE_p_uint64_t.getCPtr(fh_val), k, char_out, char_in), true);
  }

  public static SWIGTYPE_p_uint64_t ntr64(SWIGTYPE_p_uint64_t rh_val, long k, short char_out, short char_in) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.ntr64__SWIG_1(SWIGTYPE_p_uint64_t.getCPtr(rh_val), k, char_out, char_in), true);
  }

  public static SWIGTYPE_p_uint64_t ntc64(String kmer_seq, long k) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.ntc64__SWIG_0(kmer_seq, k), true);
  }

  public static SWIGTYPE_p_uint64_t ntc64(String kmer_seq, long k, SWIGTYPE_p_uint64_t fh_val, SWIGTYPE_p_uint64_t rh_val) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.ntc64__SWIG_1(kmer_seq, k, SWIGTYPE_p_uint64_t.getCPtr(fh_val), SWIGTYPE_p_uint64_t.getCPtr(rh_val)), true);
  }

  public static SWIGTYPE_p_uint64_t ntc64(short char_out, short char_in, long k, SWIGTYPE_p_uint64_t fh_val, SWIGTYPE_p_uint64_t rh_val) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.ntc64__SWIG_2(char_out, char_in, k, SWIGTYPE_p_uint64_t.getCPtr(fh_val), SWIGTYPE_p_uint64_t.getCPtr(rh_val)), true);
  }

  public static SWIGTYPE_p_uint64_t ntf64l(SWIGTYPE_p_uint64_t rh_val, long k, short char_out, short char_in) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.ntf64l(SWIGTYPE_p_uint64_t.getCPtr(rh_val), k, char_out, char_in), true);
  }

  public static SWIGTYPE_p_uint64_t ntr64l(SWIGTYPE_p_uint64_t fh_val, long k, short char_out, short char_in) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.ntr64l(SWIGTYPE_p_uint64_t.getCPtr(fh_val), k, char_out, char_in), true);
  }

  public static SWIGTYPE_p_uint64_t ntc64l(short char_out, short char_in, long k, SWIGTYPE_p_uint64_t fh_val, SWIGTYPE_p_uint64_t rh_val) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.ntc64l(char_out, char_in, k, SWIGTYPE_p_uint64_t.getCPtr(fh_val), SWIGTYPE_p_uint64_t.getCPtr(rh_val)), true);
  }

  public static SWIGTYPE_p_uint64_t ntf64(String kmer_seq, long k, long seed) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.ntf64__SWIG_2(kmer_seq, k, seed), true);
  }

  public static SWIGTYPE_p_uint64_t ntc64(String kmer_seq, long k, long seed) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.ntc64__SWIG_3(kmer_seq, k, seed), true);
  }

  public static void ntm64(String kmer_seq, long k, long m, SWIGTYPE_p_uint64_t h_val) {
    btllibJNI.ntm64__SWIG_0(kmer_seq, k, m, SWIGTYPE_p_uint64_t.getCPtr(h_val));
  }

  public static SWIGTYPE_p_uint64_t nte64(SWIGTYPE_p_uint64_t h_val, long k, long i) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.nte64(SWIGTYPE_p_uint64_t.getCPtr(h_val), k, i), true);
  }

  public static void ntm64(short char_out, short char_in, long k, long m, SWIGTYPE_p_uint64_t h_val) {
    btllibJNI.ntm64__SWIG_1(char_out, char_in, k, m, SWIGTYPE_p_uint64_t.getCPtr(h_val));
  }

  public static void ntmc64(String kmer_seq, long k, long m, SWIGTYPE_p_uint64_t h_val) {
    btllibJNI.ntmc64__SWIG_0(kmer_seq, k, m, SWIGTYPE_p_uint64_t.getCPtr(h_val));
  }

  public static void ntmc64(String kmer_seq, long k, long m, SWIGTYPE_p_uint64_t fh_val, SWIGTYPE_p_uint64_t rh_val, SWIGTYPE_p_uint64_t h_val) {
    btllibJNI.ntmc64__SWIG_1(kmer_seq, k, m, SWIGTYPE_p_uint64_t.getCPtr(fh_val), SWIGTYPE_p_uint64_t.getCPtr(rh_val), SWIGTYPE_p_uint64_t.getCPtr(h_val));
  }

  public static void ntmc64(short char_out, short char_in, long k, long m, SWIGTYPE_p_uint64_t fh_val, SWIGTYPE_p_uint64_t rh_val, SWIGTYPE_p_uint64_t h_val) {
    btllibJNI.ntmc64__SWIG_2(char_out, char_in, k, m, SWIGTYPE_p_uint64_t.getCPtr(fh_val), SWIGTYPE_p_uint64_t.getCPtr(rh_val), SWIGTYPE_p_uint64_t.getCPtr(h_val));
  }

  public static boolean ntc64(String kmer_seq, long k, SWIGTYPE_p_uint64_t h_val, SWIGTYPE_p_unsigned_int loc_n) {
    return btllibJNI.ntc64__SWIG_4(kmer_seq, k, SWIGTYPE_p_uint64_t.getCPtr(h_val), SWIGTYPE_p_unsigned_int.getCPtr(loc_n));
  }

  public static boolean ntmc64(String kmer_seq, long k, long m, SWIGTYPE_p_unsigned_int loc_n, SWIGTYPE_p_uint64_t h_val) {
    return btllibJNI.ntmc64__SWIG_3(kmer_seq, k, m, SWIGTYPE_p_unsigned_int.getCPtr(loc_n), SWIGTYPE_p_uint64_t.getCPtr(h_val));
  }

  public static boolean ntc64(String kmer_seq, long k, SWIGTYPE_p_uint64_t fh_val, SWIGTYPE_p_uint64_t rh_val, SWIGTYPE_p_uint64_t h_val, SWIGTYPE_p_unsigned_int loc_n) {
    return btllibJNI.ntc64__SWIG_5(kmer_seq, k, SWIGTYPE_p_uint64_t.getCPtr(fh_val), SWIGTYPE_p_uint64_t.getCPtr(rh_val), SWIGTYPE_p_uint64_t.getCPtr(h_val), SWIGTYPE_p_unsigned_int.getCPtr(loc_n));
  }

  public static boolean ntmc64(String kmer_seq, long k, long m, SWIGTYPE_p_uint64_t fh_val, SWIGTYPE_p_uint64_t rh_val, SWIGTYPE_p_unsigned_int loc_n, SWIGTYPE_p_uint64_t h_val) {
    return btllibJNI.ntmc64__SWIG_4(kmer_seq, k, m, SWIGTYPE_p_uint64_t.getCPtr(fh_val), SWIGTYPE_p_uint64_t.getCPtr(rh_val), SWIGTYPE_p_unsigned_int.getCPtr(loc_n), SWIGTYPE_p_uint64_t.getCPtr(h_val));
  }

  public static boolean ntmc64(String kmer_seq, long k, long m, SWIGTYPE_p_uint64_t fh_val, SWIGTYPE_p_uint64_t rh_val, SWIGTYPE_p_unsigned_int loc_n, SWIGTYPE_p_uint64_t h_val, SWIGTYPE_p_bool h_stn) {
    return btllibJNI.ntmc64__SWIG_5(kmer_seq, k, m, SWIGTYPE_p_uint64_t.getCPtr(fh_val), SWIGTYPE_p_uint64_t.getCPtr(rh_val), SWIGTYPE_p_unsigned_int.getCPtr(loc_n), SWIGTYPE_p_uint64_t.getCPtr(h_val), SWIGTYPE_p_bool.getCPtr(h_stn));
  }

  public static void ntmc64(short char_out, short char_in, long k, long m, SWIGTYPE_p_uint64_t fh_val, SWIGTYPE_p_uint64_t rh_val, SWIGTYPE_p_uint64_t h_val, SWIGTYPE_p_bool h_stn) {
    btllibJNI.ntmc64__SWIG_6(char_out, char_in, k, m, SWIGTYPE_p_uint64_t.getCPtr(fh_val), SWIGTYPE_p_uint64_t.getCPtr(rh_val), SWIGTYPE_p_uint64_t.getCPtr(h_val), SWIGTYPE_p_bool.getCPtr(h_stn));
  }

  public static SWIGTYPE_p_uint64_t mask_hash(SWIGTYPE_p_uint64_t fk_val, SWIGTYPE_p_uint64_t rk_val, String seed_seq, String kmer_seq, long k) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.mask_hash(SWIGTYPE_p_uint64_t.getCPtr(fk_val), SWIGTYPE_p_uint64_t.getCPtr(rk_val), seed_seq, kmer_seq, k), true);
  }

  public static SWIGTYPE_p_uint64_t nts64(String kmer_seq, SWIGTYPE_p_std__vectorT_bool_t seed, long k, SWIGTYPE_p_uint64_t h_val) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.nts64__SWIG_0(kmer_seq, SWIGTYPE_p_std__vectorT_bool_t.getCPtr(seed), k, SWIGTYPE_p_uint64_t.getCPtr(h_val)), true);
  }

  public static SWIGTYPE_p_uint64_t nts64(String kmer_seq, SWIGTYPE_p_std__vectorT_bool_t seed, short char_out, short char_in, long k, SWIGTYPE_p_uint64_t h_val) {
    return new SWIGTYPE_p_uint64_t(btllibJNI.nts64__SWIG_1(kmer_seq, SWIGTYPE_p_std__vectorT_bool_t.getCPtr(seed), char_out, char_in, k, SWIGTYPE_p_uint64_t.getCPtr(h_val)), true);
  }

  public static boolean ntms64(String kmer_seq, SWIGTYPE_p_std__vectorT_std__vectorT_unsigned_int_t_t seed_seq, long k, long m, SWIGTYPE_p_uint64_t fh_val, SWIGTYPE_p_uint64_t rh_val, SWIGTYPE_p_unsigned_int loc_n, SWIGTYPE_p_uint64_t h_val, SWIGTYPE_p_bool h_stn) {
    return btllibJNI.ntms64__SWIG_0(kmer_seq, SWIGTYPE_p_std__vectorT_std__vectorT_unsigned_int_t_t.getCPtr(seed_seq), k, m, SWIGTYPE_p_uint64_t.getCPtr(fh_val), SWIGTYPE_p_uint64_t.getCPtr(rh_val), SWIGTYPE_p_unsigned_int.getCPtr(loc_n), SWIGTYPE_p_uint64_t.getCPtr(h_val), SWIGTYPE_p_bool.getCPtr(h_stn));
  }

  public static void ntms64(String kmer_seq, SWIGTYPE_p_std__vectorT_std__vectorT_unsigned_int_t_t seed_seq, short char_out, short char_in, long k, long m, SWIGTYPE_p_uint64_t fh_val, SWIGTYPE_p_uint64_t rh_val, SWIGTYPE_p_uint64_t h_val, SWIGTYPE_p_bool h_stn) {
    btllibJNI.ntms64__SWIG_1(kmer_seq, SWIGTYPE_p_std__vectorT_std__vectorT_unsigned_int_t_t.getCPtr(seed_seq), char_out, char_in, k, m, SWIGTYPE_p_uint64_t.getCPtr(fh_val), SWIGTYPE_p_uint64_t.getCPtr(rh_val), SWIGTYPE_p_uint64_t.getCPtr(h_val), SWIGTYPE_p_bool.getCPtr(h_stn));
  }

  public static boolean ntmsm64(String kmer_seq, SWIGTYPE_p_std__vectorT_std__vectorT_unsigned_int_t_t seed_seq, long k, long m, long m2, SWIGTYPE_p_uint64_t fh_val, SWIGTYPE_p_uint64_t rh_val, SWIGTYPE_p_unsigned_int loc_n, SWIGTYPE_p_uint64_t h_val) {
    return btllibJNI.ntmsm64__SWIG_0(kmer_seq, SWIGTYPE_p_std__vectorT_std__vectorT_unsigned_int_t_t.getCPtr(seed_seq), k, m, m2, SWIGTYPE_p_uint64_t.getCPtr(fh_val), SWIGTYPE_p_uint64_t.getCPtr(rh_val), SWIGTYPE_p_unsigned_int.getCPtr(loc_n), SWIGTYPE_p_uint64_t.getCPtr(h_val));
  }

  public static void ntmsm64(String kmer_seq, SWIGTYPE_p_std__vectorT_std__vectorT_unsigned_int_t_t seed_seq, short char_out, short char_in, long k, long m, long m2, SWIGTYPE_p_uint64_t fh_val, SWIGTYPE_p_uint64_t rh_val, SWIGTYPE_p_uint64_t h_val) {
    btllibJNI.ntmsm64__SWIG_1(kmer_seq, SWIGTYPE_p_std__vectorT_std__vectorT_unsigned_int_t_t.getCPtr(seed_seq), char_out, char_in, k, m, m2, SWIGTYPE_p_uint64_t.getCPtr(fh_val), SWIGTYPE_p_uint64_t.getCPtr(rh_val), SWIGTYPE_p_uint64_t.getCPtr(h_val));
  }

  public static SWIGTYPE_p_std__vectorT_std__vectorT_unsigned_int_t_t parse_seeds(SWIGTYPE_p_std__vectorT_std__string_t seed_strings) {
    return new SWIGTYPE_p_std__vectorT_std__vectorT_unsigned_int_t_t(btllibJNI.parse_seeds(SWIGTYPE_p_std__vectorT_std__string_t.getCPtr(seed_strings)), true);
  }

  public static int getPIPE_READ_END() {
    return btllibJNI.PIPE_READ_END_get();
  }

  public static int getPIPE_WRITE_END() {
    return btllibJNI.PIPE_WRITE_END_get();
  }

  public static int getCOMM_BUFFER_SIZE() {
    return btllibJNI.COMM_BUFFER_SIZE_get();
  }

  public static SWIGTYPE_p_mode_t getPIPE_PERMISSIONS() {
    return new SWIGTYPE_p_mode_t(btllibJNI.PIPE_PERMISSIONS_get(), true);
  }

  public static SWIGTYPE_p_bool process_spawner_initialized() {
    return new SWIGTYPE_p_bool(btllibJNI.process_spawner_initialized(), false);
  }

  public static SWIGTYPE_p_int process_spawner_parent2child_fd() {
    long cPtr = btllibJNI.process_spawner_parent2child_fd();
    return (cPtr == 0) ? null : new SWIGTYPE_p_int(cPtr, false);
  }

  public static SWIGTYPE_p_int process_spawner_child2parent_fd() {
    long cPtr = btllibJNI.process_spawner_child2parent_fd();
    return (cPtr == 0) ? null : new SWIGTYPE_p_int(cPtr, false);
  }

  public static SWIGTYPE_p_std__mutex process_spawner_comm_mutex() {
    return new SWIGTYPE_p_std__mutex(btllibJNI.process_spawner_comm_mutex(), false);
  }

  public static long new_pipe_id() {
    return btllibJNI.new_pipe_id();
  }

  public static SWIGTYPE_p_std__mapT_std__string_btllib___Pipeline_t pipeline_map() {
    return new SWIGTYPE_p_std__mapT_std__string_btllib___Pipeline_t(btllibJNI.pipeline_map(), false);
  }

  public static String get_pipepath(long id) {
    return btllibJNI.get_pipepath(id);
  }

  public static void read_from_child(SWIGTYPE_p_void buf, long count) {
    btllibJNI.read_from_child(SWIGTYPE_p_void.getCPtr(buf), count);
  }

  public static void write_to_child(SWIGTYPE_p_void buf, long count) {
    btllibJNI.write_to_child(SWIGTYPE_p_void.getCPtr(buf), count);
  }

  public static void check_children_failures() {
    btllibJNI.check_children_failures();
  }

  public static void end_child() {
    btllibJNI.end_child();
  }

  public static void read_from_parent(SWIGTYPE_p_void buf, long count) {
    btllibJNI.read_from_parent(SWIGTYPE_p_void.getCPtr(buf), count);
  }

  public static void write_to_parent(SWIGTYPE_p_void buf, long count) {
    btllibJNI.write_to_parent(SWIGTYPE_p_void.getCPtr(buf), count);
  }

  public static boolean process_spawner_init() {
    return btllibJNI.process_spawner_init();
  }

  public static boolean getPROCESS_SPAWNER_INITIALIZER() {
    return btllibJNI.PROCESS_SPAWNER_INITIALIZER_get();
  }

  public static String get_pipeline_cmd(String path, DataStream.Operation op) {
    return btllibJNI.get_pipeline_cmd(path, op.swigValue());
  }

  public static _Pipeline run_pipeline_cmd(String cmd, DataStream.Operation op, int pipe_fd) {
    return new _Pipeline(btllibJNI.run_pipeline_cmd(cmd, op.swigValue(), pipe_fd), true);
  }

  public static SWIGTYPE_p_std__vectorT_std__string_t split(String s, String delim) {
    return new SWIGTYPE_p_std__vectorT_std__string_t(btllibJNI.split(s, delim), true);
  }

  public static void ltrim(SWIGTYPE_p_std__string s) {
    btllibJNI.ltrim(SWIGTYPE_p_std__string.getCPtr(s));
  }

  public static void rtrim(SWIGTYPE_p_std__string s) {
    btllibJNI.rtrim(SWIGTYPE_p_std__string.getCPtr(s));
  }

  public static void trim(SWIGTYPE_p_std__string s) {
    btllibJNI.trim(SWIGTYPE_p_std__string.getCPtr(s));
  }

  public static boolean starts_with(String s, String prefix) {
    return btllibJNI.starts_with(s, prefix);
  }

  public static boolean ends_with(String s, String suffix) {
    return btllibJNI.ends_with(s, suffix);
  }

  public static String get_time() {
    return btllibJNI.get_time();
  }

  public static void log_info(String msg) {
    btllibJNI.log_info(msg);
  }

  public static void log_warning(String msg) {
    btllibJNI.log_warning(msg);
  }

  public static void log_error(String msg) {
    btllibJNI.log_error(msg);
  }

  public static void check_error(boolean condition, String msg) {
    btllibJNI.check_error(condition, msg);
  }

  public static void check_warning(boolean condition, String msg) {
    btllibJNI.check_warning(condition, msg);
  }

  public static void check_stream(SWIGTYPE_p_std__ios stream, String name) {
    btllibJNI.check_stream(SWIGTYPE_p_std__ios.getCPtr(stream), name);
  }

  public static String getCOMPLEMENTS() {
    return btllibJNI.COMPLEMENTS_get();
  }

  public static String getCAPITALS() {
    return btllibJNI.CAPITALS_get();
  }

  public static void reverse_complement(SWIGTYPE_p_std__string seq) {
    btllibJNI.reverse_complement(SWIGTYPE_p_std__string.getCPtr(seq));
  }

  public static String get_reverse_complement(String seq) {
    return btllibJNI.get_reverse_complement(seq);
  }

  public static SWIGTYPE_p_unsigned_char getBIT_MASKS() {
    long cPtr = btllibJNI.BIT_MASKS_get();
    return (cPtr == 0) ? null : new SWIGTYPE_p_unsigned_char(cPtr, false);
  }

  public static long pop_cnt_byte(short x) {
    return btllibJNI.pop_cnt_byte(x);
  }

}
