package com.bcgsc.libbtl;
 
public class NtSequence {

    static {
        System.loadLibrary("btl_java");
    }

    public native void hello();
}