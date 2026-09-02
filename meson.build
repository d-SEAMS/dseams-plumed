project('dseams-plumed', 'cpp',
  version: '0.1.0',
  license: 'MIT',
  default_options: ['cpp_std=c++20', 'buildtype=release', 'warning_level=2'])

# The engine: an installed seams-core (pkg-config) or the wrap
yds_dep = dependency('seams-core', required: false)
if not yds_dep.found()
  yds_dep = dependency('seams-core', fallback: ['seams-core', 'yds_dep'])
endif

# PLUMED headers come from the installed plumed; the kernel resolves the
# symbols when LOAD dlopens the module, so nothing is linked here
plumed = find_program('plumed', required: true)
# An MPI-enabled plumed can print launcher notices before the answer; the
# path is the last line
plumed_inc = run_command(plumed, 'info', '--include-dir', check: true).stdout().strip().split('\n')[-1].strip()
plumed_root = run_command(plumed, 'info', '--root', check: true).stdout().strip().split('\n')[-1].strip()
message('plumed include dir: ' + plumed_inc)

# -isystem strings rather than include_directories(): the pixi environment
# that holds plumed can live inside this source tree, which meson refuses
# to reference through include_directories
mod = shared_module('dseams_plumed', 'src/DseamsCages.cpp',
  dependencies: [yds_dep],
  cpp_args: ['-isystem', plumed_inc, '-isystem', plumed_inc / 'plumed', '-D__PLUMED_HAS_DLOPEN'],
  name_prefix: 'lib',
  install: true,
  install_dir: get_option('libdir'))

pymod = import('python').find_installation('python3', required: false)
if pymod.found()
  test('driver_nucleation',
    pymod,
    args: [meson.project_source_root() / 'tests' / 'driver_check.py',
           plumed.full_path(), mod.full_path()],
    env: {'DSEAMS_PLUMED_DATA': meson.project_source_root() / 'tests'},
    timeout: 1800)
endif
