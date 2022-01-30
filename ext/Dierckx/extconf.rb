require 'mkmf'

rpath = File.expand_path( File.join(File.dirname(__FILE__),'DIERCKX') )
$LDFLAGS << " -L./DIERCKX -Wl,-rpath,#{rpath} "
# RbConfig::MAKEFILE_CONFIG['LD_RUN_PATH'] = './'
# create_makefile('Tst')

have_library('DIERCKX','curfit_')
have_library('DIERCKX','percur_')
have_library('DIERCKX','splev_')
have_library('DIERCKX','splder_')

create_makefile('Dierckx')
