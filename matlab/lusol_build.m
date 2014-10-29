function lusol_build
  %LUSOL_BUILD  generate thunk and prototype files for lusol interface

  % unload the library if it is loaded
  if libisloaded('libclusol')
    unloadlibrary('libclusol');
  end
  
  % load the library to generate proto and thunk files
  [notfound, warnings] = loadlibrary('libclusol','clusol.h','notempdir');

  if ~isempty(notfound)
    error('lusol_build:loadlibrary','functions missing from library');
  end
  
  if ~isempty(warnings)
    warning('lusol_build:loadlibrary','loadlibrary returned warning messages:\n%s\n',warnings);
  end
  
  % rename the thunk file
  movefile('libclusol_proto.m',['libclusol_proto_' lower(computer) '.m'])

  % clean up the code files
  tcfilename = ['libclusol_thunk_' lower(computer) '.c'];
  delete('clusol.i',tcfilename);

  % finally unload the library
  unloadlibrary('libclusol');

end
