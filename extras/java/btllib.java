/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.1
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package btllib;

public class btllib {
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

  public static SWIGTYPE_p_FILE data_load(String source) {
    long cPtr = btllibJNI.data_load(source);
    return (cPtr == 0) ? null : new SWIGTYPE_p_FILE(cPtr, false);
  }

  public static SWIGTYPE_p_FILE data_save(String sink, boolean append) {
    long cPtr = btllibJNI.data_save(sink, append);
    return (cPtr == 0) ? null : new SWIGTYPE_p_FILE(cPtr, false);
  }

  public static void sigchld_handler(int sig) {
    btllibJNI.sigchld_handler(sig);
  }

  public static boolean data_saveload_init() {
    return btllibJNI.data_saveload_init();
  }

  public static boolean getData_saveload_initialized() {
    return btllibJNI.data_saveload_initialized_get();
  }

  public static String get_saveload_cmd(String path, SaveloadOp op) {
    return btllibJNI.get_saveload_cmd(path, op.swigValue());
  }

  public static SWIGTYPE_p_FILE run_saveload_cmd(String cmd, SaveloadOp op) {
    long cPtr = btllibJNI.run_saveload_cmd(cmd, op.swigValue());
    return (cPtr == 0) ? null : new SWIGTYPE_p_FILE(cPtr, false);
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

}
