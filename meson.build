project('mesas', 'fortran')

python = import('python').find_installation(pure: false)

sas_solve = python.extension_module(
  'mesas.sas.solve',
  'mesas/sas/solve.f90',
  install: true,
  install_dir: python.get_install_dir(),
  fortran_args: ['-fPIC', '-fno-stack-arrays']
)

install_subdir('mesas', install_dir: python.get_install_dir())
