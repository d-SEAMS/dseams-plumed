project('dseams-plumed', 'cpp',
  version: '2.14.0',
  license: 'MIT',
  default_options: ['cpp_std=c++20', 'buildtype=release', 'warning_level=2'])

# The engine: an installed seams-core (pkg-config) or the wrap.
# The wrap pins the seams-core 2.14.0 release, the first engine release
# that distinguishes the published Nguyen and Molinero rule from
# ovito-chill+. --force-fallback-for=seams-core selects it over an older
# installed copy. The engine tests are not built. The CLI stays on because
# driver_library writes its key library with it. The engine's own LAMMPS,
# GROMACS, OpenMM and PLUMED hooks are auto features that this module
# does not use; left on, the LAMMPS hook builds whenever liblammps is
# visible and needs an MPI-enabled LAMMPS.
yds_dep = dependency('seams-core',
  fallback: ['seams-core', 'yds_dep'],
  default_options: ['with_tests=false', 'with_cli=true',
                    'with_lammps=disabled', 'with_gromacs=disabled',
                    'with_openmm=disabled', 'with_plumed=disabled'])

# PLUMED headers come from the installed plumed; libplumedKernel is linked
# below
plumed = find_program('plumed', required: true)
# An MPI-enabled plumed can print launcher notices before the answer; the
# path is the last line
plumed_inc = run_command(plumed, 'info', '--include-dir', check: true).stdout().strip().split('\n')[-1].strip()
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
# One set of flags for every target that compiles smooth_chill.hpp, so
# thread_bits and track check the code the module runs. -O3 is repeated
# because conda's CXXFLAGS end in -O2, after the buildtype's -O3.
kernel_args = ['-O3'] + (get_option('native') ? ['-march=native'] : [])
mod = shared_module('dseams_plumed',
  'plumed/dseams/DseamsCages.cpp',
  'plumed/dseams/DseamsChill.cpp',
  'plumed/dseams/DseamsRings.cpp',
  'plumed/dseams/DseamsSmooth.cpp',
  dependencies: [yds_dep, plumed_kernel, omp],
  cpp_args: ['-isystem', plumed_inc, '-isystem', plumed_inc / 'plumed', '-D__PLUMED_HAS_DLOPEN'] + kernel_args,
  name_prefix: 'lib',
  install: true,
  install_dir: get_option('libdir'))

offline = executable('dseams_offline', 'tests/offline_frame.cpp',
  dependencies: [yds_dep, omp],
  include_directories: include_directories('plumed/dseams'),
  cpp_args: kernel_args,
  install: false)

thread_bits = executable('dseams_thread_bits', 'tests/thread_bits.cpp',
  dependencies: [omp],
  include_directories: include_directories('plumed/dseams'),
  cpp_args: kernel_args,
  install: false)
test('thread_bits', thread_bits, timeout: 120)
track = executable('dseams_track', 'tests/track_check.cpp',
  dependencies: [yds_dep, omp],
  include_directories: include_directories('plumed/dseams'),
  cpp_args: kernel_args,
  install: false)
test('track', track, timeout: 180)

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
  # The library check needs the engine CLI to write the keys: the CLI built
  # by the seams-core subproject when the wrap supplies the engine, else a
  # seams on PATH. Without either the test is registered and reports a
  # skip (exit 77) with the reason.
  # The wrap is configured with with_cli=true, so seams_exe exists there.
  seams_cli = ''
  if yds_dep.type_name() == 'internal'
    seams_cli = subproject('seams-core').get_variable('seams_exe')
  else
    seams_path = find_program('seams', required: false)
    if seams_path.found()
      seams_cli = seams_path
    endif
  endif
  test('driver_library',
    pymod,
    args: [meson.project_source_root() / 'tests' / 'library_check.py',
           plumed.full_path(), mod.full_path(), seams_cli],
    timeout: 600)
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
