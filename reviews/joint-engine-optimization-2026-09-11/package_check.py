"""Stage distribution verification without overwriting existing build outputs."""
import os
import shutil
import subprocess
import sys
from pathlib import Path

source = Path('/work')
target = Path('/tmp/nwkit-package-check')
target.mkdir()
for path in source.iterdir():
    if path.is_file() and not path.name.startswith('.'):
        shutil.copy2(path, target / path.name)
ignore = shutil.ignore_patterns('__pycache__', '*.pyc')
for name in ('nwkit', 'tests', 'tools', 'examples', 'scripts'):
    if (source / name).exists():
        shutil.copytree(source / name, target / name, ignore=ignore)
for line in (source / 'MANIFEST.in').read_text().splitlines():
    parts = line.split()
    if len(parts) == 2 and parts[1].startswith('reviews/'):
        path = source / parts[1]
        output = target / parts[1]
        output.parent.mkdir(parents=True, exist_ok=True)
        if path.is_dir():
            shutil.copytree(path, output, ignore=ignore)
        elif path.exists():
            shutil.copy2(path, output)
environment = os.environ.copy()
environment['SOURCE_DATE_EPOCH'] = '1789096833'
environment['PATH'] = str(Path(sys.executable).parent) + os.pathsep + environment['PATH']
subprocess.run([sys.executable, 'tools/check.py', 'dist'], cwd=target, env=environment, check=True)

# Exercise the built source distribution and installed wheel, not this checkout.
import tarfile

extracted = Path('/tmp/nwkit-sdist-check')
extracted.mkdir()
with tarfile.open(next((target / 'dist').glob('*.tar.gz'))) as archive:
    archive.extractall(extracted, filter='data')
sdist = next(extracted.iterdir())
environment.pop('PYTHONPATH', None)
subprocess.run([sys.executable, '-m', 'pytest', 'tests/test_cli.py', 'tests/test_cli_contracts.py', 'tests/test_wiki_examples.py', '-q'], cwd=sdist, env=environment, check=True)
wheel = next((target / 'dist').glob('*.whl'))
subprocess.run([sys.executable, '-m', 'pip', 'install', '--force-reinstall', '--no-deps', str(wheel)], env=environment, check=True)
subprocess.run([str(Path(sys.executable).parent / 'nwkit'), '--version'], cwd='/tmp', env=environment, check=True)
help_result = subprocess.run([str(Path(sys.executable).parent / 'nwkit'), 'shift', '--help'], cwd='/tmp', env=environment, check=True, capture_output=True, text=True)
assert 'bounded dense GLS' in help_result.stdout
smoke = '''
import numpy as np
from ete4 import Tree
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout
import nwkit.shift_joint_dense as dense
print('Installed engine:', dense.__file__)
tree=Tree('(((a:0.4,b:0.4):0.3,(c:0.4,d:0.4):0.3):0.3,((e:0.4,f:0.4):0.3,(g:0.4,h:0.4):0.3):0.3);')
values=np.random.default_rng(71).normal(size=(8,2))
values[1,0]=values[3,1]=np.nan
data=ShiftData.build(tree,values,('x','y'),np.full((8,2),.04))
fit=fit_native_layout(data,ShiftLayout.build(data.tree),options=NativeFitOptions(trait_covariance='full',optimizer_starts=2),alpha_height=[.5,2.])
assert fit['joint_covariance']['engine']=='dense_observed_gls'
assert np.isfinite(fit['log_likelihood'])
print('Installed wheel joint GLS smoke passed')
'''
subprocess.run([sys.executable, '-c', smoke], cwd='/tmp', env=environment, check=True)
