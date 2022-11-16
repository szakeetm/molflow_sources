Works on mac and linux.
Check out repo: git clone https://bitbucket.org/wbenny/ziplib.git
If not available use zipped source
Copy modified makefiles (linux: just make sure it builds static lib, mac: remove -fPIC option)
make -f MakeFile_linux or similar
Copy resulting libzip.a to lib_external/*/libzip_gcc.a
