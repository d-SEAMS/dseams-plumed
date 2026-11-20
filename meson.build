project('dseams-plumed', 'cpp',
  version: '0.4.1',
  license: 'MIT',
  default_options: ['cpp_std=c++20', 'buildtype=release', 'warning_level=2'])

# The engine: an installed seams-core (pkg-config seams-core) at 2.14.0 or
# newer, or the wrap. The wrap revision is the tag v2.14.0.
yds_dep = dependency('seams-core', version: '>=2.14.0', required: false)
seams_from_wrap = not yds_dep.found()
if seams_from_wrap
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
    # `seams read` of a missing or unreadable PDB, GRO, DCD, or CON file
    # exits 2 with `cannot read <path>` only when those readers are in the
    # engine. Without them the CLI throws "not in this build" instead, so
    # these four tests are registered only for a wrap build that has both.
    # A build with the readers disabled still passes every test it
    # registers (the four drivers and feature_query). CI enables both.
    # An installed engine is not the patched CLI this test drives.
    engine_readers = false
    if seams_from_wrap
      engine_readers = subproject('seams-core').get_variable('seams_has_structure_readers')
    endif
    if engine_readers
      foreach ext : ['pdb', 'gro', 'dcd', 'con']
        test('reader_fail_' + ext,
          pymod,
          args: [meson.project_source_root() / 'tests' / 'reader_fail_check.py',
                 seams_cli.full_path(), ext],
          timeout: 60)
      endforeach
    endif
  endif
endif
