function lusol_build
  %LUSOL_BUILD  generate thunk and prototype files for lusol interface

  % unload the library if it is loaded
  if libisloaded('libclusol')
    unloadlibrary('libclusol');
  end

  % get the compiler setting
  cc = mex.getCompilerConfigurations('C','Selected');
  cc = cc(1);

  % copy fields of compiler configuration
  cc_field_cell = fieldnames(cc);
  for cc_field = cc_field_cell(:)'
    field_str = cc_field{1};
    cc1.(field_str) = cc.(field_str);
  end

  % remove details field
  cc1 = rmfield(cc1,'Details');

  % copy details
  cc_detail_cell = fieldnames(cc.Details);
  for cc_detail_field = cc_detail_cell(:)'
    field_str = cc_detail_field{1};
    cc1_detail.(field_str) = cc.Details.(field_str);
  end

  % set the compiler
  cc1_detail.CompilerExecutable = 'gcc';

  % create new compiler configuration
  cc = mex.CompilerConfiguration(cc1,cc1_detail);

  % load the library to generate proto and thunk files
  loadlibrary('libclusol','clusol.h','notempdir','compilerconfiguration',cc);

  % rename the thunk file
  movefile('libclusol_proto.m',['libclusol_proto_' lower(computer) '.m'])

  % clean up the code files
  tcfilename = ['libclusol_thunk_' lower(computer) '.c'];
  delete('clusol.i',tcfilename);

  % finally unload the library
  unloadlibrary('libclusol');

end
