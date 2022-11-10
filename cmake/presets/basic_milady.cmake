# preset that turns on just a few, frequently used packages
# this will be compiled quickly and handle a lot of common inputs.

set(ALL_PACKAGES EXTRA-FIX KSPACE MANYBODY MISC ML-MILADY ML-SNAP MOLECULE REPLICA RIGID)

foreach(PKG ${ALL_PACKAGES})
  message(STATUS "!!!! preset on - ${PKG}")
  set(PKG_${PKG} ON CACHE BOOL "on" FORCE)
endforeach()
