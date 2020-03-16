$nrnpath = Join-Path -Path $env:neuronhome -ChildPath "bin\nrniv.exe"
& $nrnpath -NSTACK 100000 -NFRAME 20000 -Py_NoSiteFlag init_field_stim_test.hoc