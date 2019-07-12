package com.bcgsc.libbtl;
 
public class Sequence {

    static {
        System.loadLibrary("btl_java");
    }

    public native void hello();
}