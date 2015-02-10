require 'mkmf'

rpath = File.expand_path( File.join(File.dirname(__FILE__), 'LIMEX4_3A') )
$LDFLAGS << " -L./LIMEX4_3A -Wl,-rpath,#{rpath} "
# RbConfig::MAKEFILE_CONFIG['LD_RUN_PATH'] = './'
# create_makefile('Tst')

have_library('LIMEX4_3A', 'limex_')
# have_library('ODEydot', 'ydot_limex_')
# have_library('ODEydot', 'init_ode_') 
# #have_library('ODEydot', 'get_species_ids_')
# #have_library('ODEydot', 'get_parameter_ids_')

create_makefile('Limex')
# create_makefile('LimexDL')
