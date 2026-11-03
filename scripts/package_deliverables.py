"""Assemble self-contained, credential-free project downloads."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from shutil import copy2
ROOT=Path(__file__).resolve().parents[1]
assets=ROOT/'website/dist/assets'
for name in ['report.docx','report.pdf']:
 copy2(ROOT/'deliverables'/name,assets/name)
for name in ['results.csv','diagnostics.csv','summary.csv','paired_differences.csv','manifest.json']:
 copy2(ROOT/'outputs/main'/name,assets/name)
files=[ROOT/name for name in ['model.py','run_experiments.py','figures.py','supplement.py','supplement_figures.py','README.md','RESULTS.md','requirements.txt','requirements-lock.txt','frozen_config.json']]
files+=list((ROOT/'tests').glob('*.py'))+list((ROOT/'scripts').glob('*.py'))
files+=[ROOT/'deliverables/report.md',ROOT/'SUBMISSION_CHECKLIST.md',ROOT/'236 proposal.pdf']
files+=list((ROOT/'data/example').glob('*'))
files+=list((ROOT/'.github/workflows').glob('*.yml'))
for stage in ['pilot','main']:
 files+=list((ROOT/'outputs'/stage).glob('*.csv'))+list((ROOT/'outputs'/stage).glob('*.json'))
 files+=list((ROOT/'outputs'/stage/'datasets').glob('*.json'))
 files+=list((ROOT/'outputs'/stage/'figures').glob('*.png'))
files+=list((ROOT/'report-assets').glob('*'))
files += [p for p in (ROOT/'outputs/supplement').rglob('*') if p.is_file() and p.suffix in {'.csv','.json','.png'}]
files += [p for p in (assets/'supplement').rglob('*') if p.is_file()]
# Include portable website source, excluding binary downloads to avoid recursion.
files+=list((ROOT/'website/dist').glob('*.html'))+list((ROOT/'website/dist').glob('*.css'))+list((ROOT/'website/dist').glob('*.js'))
files+=[assets/'results.json',assets/'network_example.png',assets/'favicon.svg']
archive=ROOT/'deliverables/reproducibility.zip'
with ZipFile(archive,'w',ZIP_DEFLATED) as z:
 for path in sorted(set(files)):
  z.write(path,Path('stat236-project')/path.relative_to(ROOT))
copy2(archive,assets/archive.name)
with ZipFile(ROOT/'deliverables/website.zip','w',ZIP_DEFLATED) as z:
 for path in sorted((ROOT/'website/dist').rglob('*')):
  if path.is_file():z.write(path,path.relative_to(ROOT/'website/dist'))
print('Created reproducibility.zip and website.zip without virtual environments or hosting credentials.')
