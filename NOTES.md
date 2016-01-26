# LUSOL-Matlab interface notes

This document contains links to references on building shared/dynamic libraries on Mac and Linux.

## Mac OS X references

Apple Developer Documentation: <https://developer.apple.com/library/mac/navigation/>

Specific documents:

* [Dynamic Library Programming Topics][dlpt]
* [Porting UNIX/Linux Applications to OS X][pula]

[dlpt]: https://developer.apple.com/library/mac/documentation/DeveloperTools/Conceptual/DynamicLibraries/000-Introduction/Introduction.html#//apple_ref/doc/uid/TP40001908-SW1
[pula]: https://developer.apple.com/library/mac/documentation/Porting/Conceptual/PortingUnix/intro/intro.html#//apple_ref/doc/uid/TP40002847-TPXREF101

## Linux references

### **Linux Programming Interface**:

* <http://proquest.safaribooksonline.com/9781593272203?uicode=stanford>
* See chapters on shared libraries (41 and 42)

### TLDP: Program Library HOWTO

<http://tldp.org/HOWTO/Program-Library-HOWTO/index.html>

### How to write shared libraries by Drepper

<https://www.akkadia.org/drepper/dsohowto.pdf>

### Tracking down libraries used

```
$ otool -L libgfortran.dylib 
libgfortran.dylib:
	/usr/local/opt/gcc/lib/gcc/5/libgfortran.3.dylib (compatibility version 4.0.0, current version 4.0.0)
	/usr/local/Cellar/gcc/5.3.0/lib/gcc/5/libquadmath.0.dylib (compatibility version 1.0.0, current version 1.0.0)
	/usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1225.1.1)
	/usr/local/lib/gcc/5/libgcc_s.1.dylib (compatibility version 1.0.0, current
    version 1.0.0)
```

```
$ otool -L libquadmath.dylib 
libquadmath.dylib:
	/usr/local/opt/gcc/lib/gcc/5/libquadmath.0.dylib (compatibility version 1.0.0, current version 1.0.0)
	/usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1225.1.1)
	/usr/local/lib/gcc/5/libgcc_s.1.dylib (compatibility version 1.0.0, current
    version 1.0.0)
```

LDFLAGS += /usr/local/opt/gcc/lib/gcc/5/libgfortran.a
LDFLAGS += /usr/local/opt/gcc/lib/gcc/5/libquadmath.a
LDFLAGS += -L/usr/local/lib/gcc/5 -lgcc_s.1
