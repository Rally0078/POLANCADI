# POLAN for IDL CADI software

POLAN with modified input format to be compatible with the modified cadi_igwscale.pro to compute real heights from virtual heights on an ionogram.

An example of the input format is shown in the textfile in.dat. The IDL script cadi_igwscale.pro included in the repo will output data of that format into the in.dat file. Default locations of the files are:

The POLAN executable
```
C:\polan\polan.exe
```

in.dat
```
C:\polan\in.dat
```

cadi_igwscale.pro
```
C:\idl\IFG_CADI\cadi_igwscale.pro
```

The output of POLAN is stored in out.dat, and the format remains unchanged from the original POLAN.