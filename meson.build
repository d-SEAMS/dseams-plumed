project('dseams-plumed', 'cpp',
  version: '0.4.1',
  license: 'MIT',
  default_options: ['cpp_std=c++20', 'buildtype=release', 'warning_level=2'])

# The engine: an installed seams-core (pkg-config seams-core) at 2.14.0 or
# newer, or the wrap. The wrap revision is the tag v2.14.0.
yds_dep = dependency('seams-core', version: '>=2.14.0', required: false)
if not yds_dep.found()
  yds_dep = dependency('seams-core', version: '>=2.14.0',
    fallback: ['seams-core', 'yds_dep'])
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
  # Reports dseams_plumed_has_chemfiles()/has_readcon() every build, whether
  # or not the engine has those optional readers: a build without them must
  # still succeed (and this test with it), only the printed values differ.
  # ci.yml's build-test job configures with -Dseams-core:with_chemfiles=enabled
  # -Dseams-core:with_readcon=enabled explicitly (a deliberate configuration,
  # not auto-detection), so it additionally runs this script with
  # --require-chemfiles --require-readcon, which fails (not skips) if either
  # is unexpectedly missing there.
  test('feature_query',
    pymod,
    args: [meson.project_source_root() / 'tests' / 'feature_check.py',
           plumed.full_path(), mod.full_path()],
    timeout: 60)
  # the library check needs the engine CLI to write the keys. An
  # installed seams is used when it is on PATH; otherwise the CLI the
  # wrap just built.
  seams_cli = find_program('seams', required: false)
  if not seams_cli.found()
    seams_cli = subproject('seams-core').get_variable('seams_exe')
  endif
  if seams_cli.found()
    test('driver_library',
      pymod,
      args: [meson.project_source_root() / 'tests' / 'library_check.py',
             plumed.full_path(), mod.full_path(), seams_cli.full_path()],
      timeout: 600)
  endif
endif
