require 'mkmf'

rpath = File.expand_path( File.join(File.dirname(__FILE__),'NLSCON') )
$LDFLAGS << " -L./NLSCON -Wl,-rpath,#{rpath} "
# RbConfig::MAKEFILE_CONFIG['LD_RUN_PATH'] = './'
# create_makefile('Tst')

have_library('NLSCON','nlscon_')

create_makefile('Nlscon')
