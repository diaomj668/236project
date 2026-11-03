"""Create the editable report using the bundled document runtime."""
from pathlib import Path
import re
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
ROOT=Path(__file__).resolve().parents[1]
source=ROOT/'deliverables/report.md'
text=source.read_text().replace('{{ALL_ERRORS}}',(ROOT/'report-assets/all_errors.md').read_text())
source.write_text(text)
doc=Document(); sec=doc.sections[0]
sec.page_width=Inches(8.5);sec.page_height=Inches(11)
sec.top_margin=Inches(.7);sec.bottom_margin=Inches(.7)
sec.left_margin=Inches(.8);sec.right_margin=Inches(.8)
for name in ['Normal','Title','Subtitle','Heading 1','Heading 2','Heading 3','Caption']:
 st=doc.styles[name];st.font.name='Times New Roman';st.font.color.rgb=RGBColor(0,0,0)
 st.paragraph_format.space_after=Pt(6);st.paragraph_format.line_spacing=1
st=doc.styles['Normal'];st.font.size=Pt(12);st.paragraph_format.widow_control=True
for name,size in [('Title',22),('Heading 1',15),('Heading 2',13),('Heading 3',12),('Caption',10)]:
 doc.styles[name].font.size=Pt(size)
 doc.styles[name].paragraph_format.space_before=Pt(10 if name!='Title' else 0)
for style in doc.styles:
 if style.type == 1:
  props=style.element.get_or_add_pPr()
  for e in list(props):
   if e.tag==qn('w:pBdr'):props.remove(e)
  rpr=style.element.get_or_add_rPr()
  fonts=rpr.find(qn('w:rFonts'))
  if fonts is not None:
   for key in list(fonts.attrib):
    if 'Theme' in key:del fonts.attrib[key]
   fonts.set(qn('w:ascii'),'Times New Roman');fonts.set(qn('w:hAnsi'),'Times New Roman')
doc.styles['Caption'].font.bold=False
header=sec.header.paragraphs[0];header.text='STAT 236  |  From Purchases to Communities'
header.runs[0].font.size=Pt(9)
footer=sec.footer.paragraphs[0];footer.alignment=WD_ALIGN_PARAGRAPH.CENTER
field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');footer._p.append(field)

def add_inline(p,text):
 for token in re.split(r'(\*\*.*?\*\*|`[^`]+`)',text):
  if token.startswith('**'):
   p.add_run(token[2:-2]).bold=True
  elif token.startswith('`'):
   run=p.add_run(token[1:-1]);run.font.name='Courier New';run.font.size=Pt(10)
  else:p.add_run(token)

def add_table(lines):
 rows=[[c.strip() for c in line.strip().strip('|').split('|')] for line in lines]
 rows=[r for r in rows if not all(re.fullmatch('[-: ]+',c) for c in r)]
 table=doc.add_table(rows=0,cols=len(rows[0]));table.alignment=WD_TABLE_ALIGNMENT.CENTER
 table.autofit=False
 widths=([.65,.9,1.4,1,1,1, .95] if len(rows[0])==7 else
         [.9,1.85,2.05,2.1] if len(rows[0])==4 and rows[0][0]=='Mean purchases' else
         [1.2,.8,2,2.9])
 if rows[0][0]=='Setting':widths=[1.,.5,.5,1.2,1.2,1.2,1.2]
 if rows[0][0]=='Cell':widths=[1.55,.85,.7,.8,.7,1.1,1.25]
 for j,col in enumerate(table.columns):col.width=Inches(widths[j])
 for i,row in enumerate(rows):
  cells=table.add_row().cells
  for j,(cell,value) in enumerate(zip(cells,row)):
   cell.width=Inches(widths[j]);cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
   p=cell.paragraphs[0];p.paragraph_format.space_after=Pt(4);p.paragraph_format.space_before=Pt(4)
   p.alignment=WD_ALIGN_PARAGRAPH.CENTER
   r=p.add_run(value);r.font.size=Pt(10 if len(row)==7 else 10.5);r.bold=i==0
   props=cell._tc.get_or_add_tcPr();borders=OxmlElement('w:tcBorders')
   for edge in ['top','left','bottom','right']:
    e=OxmlElement('w:'+edge);e.set(qn('w:val'),'single');e.set(qn('w:sz'),'4');e.set(qn('w:color'),'D9D9D9');borders.append(e)
   props.append(borders)
   sh=OxmlElement('w:shd');sh.set(qn('w:fill'),'E8EDF2' if i==0 else 'FFFFFF');props.append(sh)
  pr=table.rows[-1]._tr.get_or_add_trPr();n=OxmlElement('w:cantSplit');pr.append(n)
  if i==0:
   repeat=OxmlElement('w:tblHeader');pr.append(repeat)
 sp=doc.add_paragraph();sp.paragraph_format.space_after=Pt(0);sp.paragraph_format.space_before=Pt(0);sp.paragraph_format.line_spacing=Pt(1);sp.add_run().font.size=Pt(1)

lines=text.splitlines();i=0;in_code=False;new_page=False
while i<len(lines):
 line=lines[i].strip()
 if not line:i+=1;continue
 if line=='<!-- page -->':new_page=True;i+=1;continue
 if line.startswith('```'):in_code=not in_code;i+=1;continue
 if in_code:
  p=doc.add_paragraph();p.paragraph_format.space_after=Pt(3)
  r=p.add_run(line);r.font.name='Courier New';r.font.size=Pt(9)
  i+=1;continue
 if line.startswith('|'):
  group=[]
  while i<len(lines) and lines[i].strip().startswith('|'):group.append(lines[i]);i+=1
  add_table(group);continue
 if line.startswith('!['):
  match=re.match(r'!\[(.*?)\]\((.*?)\)',line)
  p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.paragraph_format.keep_with_next=True
  p.add_run().add_picture(str((source.parent/match[2]).resolve()),width=Inches(5.3 if 'network_example' in match[2] else 6.3 if 'sensitivity_r3' in match[2] else 6.8))
  i+=1;continue
 if line.startswith('$$'):
  p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.paragraph_format.space_after=Pt(7)
  r=p.add_run(line.strip('$ ').strip());r.italic=True;r.font.size=Pt(11)
  i+=1;continue
 if line.startswith('# '):doc.add_paragraph(line[2:],'Title');i+=1;continue
 if line.startswith('## '):doc.add_paragraph(line[3:],'Subtitle');i+=1;continue
 if line.startswith('### '):
  p=doc.add_paragraph(line[4:],'Heading 1')
  if new_page:p.paragraph_format.page_break_before=True;new_page=False
  i+=1;continue
 if re.match(r'^(Figure|Table) [0-9]+\.',line):
  doc.add_paragraph(line,'Caption');i+=1;continue
 p=doc.add_paragraph();add_inline(p,line);i+=1
for p in doc.paragraphs:
 for e in list(p._p.get_or_add_pPr()):
  if e.tag==qn('w:pBdr'):p._p.get_or_add_pPr().remove(e)
doc.core_properties.title='From Purchases to Communities'
doc.core_properties.author='Tong Zhou; Sia Xie; Yilin Jiang; Mingjiao Diao'
doc.core_properties.subject='STAT 236 simulation study of activity heterogeneity and preference breadth'
doc.save(ROOT/'deliverables/report.docx')
print(ROOT/'deliverables/report.docx')
