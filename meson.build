project('dseams-plumed', 'cpp',
  version: '0.3.1',
  license: 'MIT',
  default_options: ['cpp_std=c++20', 'buildtype=release', 'warning_level=2'])

# The engine: an installed seams-core (pkg-config) or the wrap.
# The wrap is HaoZeke/seams-core chill-plus/published-rule-pr, which is the
# tree that distinguishes the published Nguyen and Molinero rule from
# ovito-chill+. --force-fallback-for=seams-core selects it over an older
# installed copy. The engine tests are not built. The CLI stays on because
# the language bindings refer to it.
yds_dep = dependency('seams-core',
  fallback: ['seams-core', 'yds_dep'],
  default_options: ['with_tests=false', 'with_cli=true'])

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
# Link the kernel as well as including its headers. LAMMPS loads the kernel
# into the global namespace, so an undefined PLUMED symbol in the module is
# enough. GROMACS dlopens the kernel locally, and the module then has to
# name libplumedKernel itself.
plumed_kernel = meson.get_compiler('cpp').find_library('plumedKernel', required: true)
# The smooth fraction uses the same OpenMP width as PLUMED (NUM_THREADS).
omp = dependency('openmp', required: true)
mod = shared_module('dseams_plumed',
  'src/DseamsCages.cpp',
  'src/DseamsChill.cpp',
  'src/DseamsRings.cpp',
  'src/DseamsSmooth.cpp',
  dependencies: [yds_dep, plumed_kernel, omp],
  cpp_args: ['-isystem', plumed_inc, '-isystem', plumed_inc / 'plumed', '-D__PLUMED_HAS_DLOPEN'],
  name_prefix: 'lib',
  install: true,
  install_dir: get_option('libdir'))

offline = executable('dseams_offline', 'tests/offline_frame.cpp',
  dependencies: [yds_dep, omp],
  include_directories: include_directories('src'),
  install: false)

thread_bits = executable('dseams_thread_bits', 'tests/thread_bits.cpp',
  dependencies: [omp],
  include_directories: include_directories('src'),
  install: false)
test('thread_bits', thread_bits, timeout: 120)

pymod = import('python').find_installation('python3', required: false)
if pymod.found()
  test('driver_nucleation',
    pymod,
    args: [meson.project_source_root() / 'tests' / 'driver_check.py',
           plumed.full_path(), mod.full_path()],
    env: {'DSEAMS_PLUMED_DATA': meson.project_source_root() / 'tests'},
    timeout: 1800)
  test('driver_ions',
    pymod,
    args: [meson.project_source_root() / 'tests' / 'ion_check.py',
           plumed.full_path(), mod.full_path()],
    timeout: 600)
  test('driver_guests',
    pymod,
    args: [meson.project_source_root() / 'tests' / 'guest_check.py',
           plumed.full_path(), mod.full_path()],
    timeout: 600)
  # the library check needs the engine CLI on PATH to write the keys
  seams_cli = find_program('seams', required: false)
  if seams_cli.found()
    test('driver_library',
      pymod,
      args: [meson.project_source_root() / 'tests' / 'library_check.py',
             plumed.full_path(), mod.full_path(), seams_cli.full_path()],
      timeout: 600)
  endif
  test('driver_chill',
    pymod,
    args: [meson.project_source_root() / 'tests' / 'chill_match.py',
           plumed.full_path(), mod.full_path(), offline.full_path()],
    timeout: 1800)
  test('regtest',
    find_program('bash'),
    args: [meson.project_source_root() / 'regtest' / 'run'],
    env: {'DSEAMS_MODULE': mod.full_path(), 'PLUMED_BIN': plumed.full_path()},
    timeout: 600)
  test('driver_smooth',
    pymod,
    args: [meson.project_source_root() / 'tests' / 'smooth_check.py',
           plumed.full_path(), mod.full_path()],
    timeout: 1800)
endif
