set destdir="source_snapshots\%date%_%time:~0,2%%time:~3,2%"
md %destdir%
xcopy /e /exclude:snapshotexclude.txt Source_files %destdir%