"""Prospective predictions and evaluation of the two targeted well points."""
import sys,json,hashlib,math
from analyze import ROOT,read,write,value
FROZEN=ROOT/'data/stage2_frozen_models.json'
PRED=ROOT/'data/well_65_70_predictions.json'
def freeze():
    models=json.loads(FROZEN.read_text())['models']
    record={'source_model_sha256':hashlib.sha256(FROZEN.read_bytes()).hexdigest(),
      'selection':'targeted near known gauche minimum; no refit',
      'predictions':[dict(order=n,angle_deg=a,predicted_kj_mol=value(a,m['coefficients_kj_mol'])) for n,m in models.items() for a in [65,70]]}
    text=json.dumps(record,indent=2)+'\n'
    if PRED.exists() and PRED.read_text()!=text:raise ValueError('Frozen predictions changed')
    PRED.write_text(text,encoding='utf-8',newline='\n')
    print(text)
def evaluate():
    record=json.loads(PRED.read_text());assert hashlib.sha256(FROZEN.read_bytes()).hexdigest()==record['source_model_sha256']
    truth={int(r['target_angle_deg']):float(r['relative_electronic_kj_mol']) for r in read('well_65_70_audit.csv')}
    assert sorted(truth)==[65,70]
    residuals=[];metrics=[]
    for p in record['predictions']:
        e=truth[p['angle_deg']];residuals.append(dict(**p,calculated_kj_mol=e,error_kj_mol=p['predicted_kj_mol']-e))
    for n in ['3','4','5']:
        errors=[r['error_kj_mol'] for r in residuals if r['order']==n]
        metrics.append(dict(order=n,n_targeted=2,mae_kj_mol=sum(map(abs,errors))/2,rmse_kj_mol=math.sqrt(sum(e*e for e in errors)/2),max_abs_kj_mol=max(map(abs,errors))))
    write('well_65_70_residuals.csv',residuals);write('well_65_70_metrics.csv',metrics)
    plot(truth,json.loads(FROZEN.read_text())['models'])
    print(json.dumps({'residuals':residuals,'metrics':metrics},indent=2))
def plot(truth,models):
    def x(a):return 70+(a-55)*24
    def y(e):return 370-(e-3)*75
    p=['<svg xmlns="http://www.w3.org/2000/svg" width="800" height="490" viewBox="0 0 800 490">','<rect width="800" height="490" fill="white"/>','<g font-family="sans-serif" font-size="14">','<text x="70" y="25" font-size="20">Gauche well: frozen curves and two new calculations</text>']
    for e in [3,4,5,6,7]:p += [f'<path d="M70 {y(e)} H670" stroke="#ddd"/>',f'<text x="35" y="{y(e)+5}">{e}</text>']
    for a in [55,60,65,70,75,80]:p.append(f'<text x="{x(a)-8}" y="400">{a}</text>')
    for i,(n,m) in enumerate(models.items()):
        color=['#777','#176cc0','#c66a18'][i]
        pts=' '.join(f'{x(55+j/10):.2f},{y(value(55+j/10,m["coefficients_kj_mol"])):.2f}' for j in range(251))
        p += [f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2" stroke-dasharray="5 4"/>',f'<text x="{70+i*210}" y="460" fill="{color}">Order {n} (frozen)</text>']
    for a,e in truth.items():p.append(f'<rect x="{x(a)-5}" y="{y(e)-5}" width="10" height="10" fill="white" stroke="#ac204a" stroke-width="2"/>')
    p += ['<text x="70" y="50">Squares: new Gaussian values; dashed curves: predictions</text>','<text x="70" y="430">Angle (degrees); relative electronic energy (kJ/mol)</text>','</g></svg>']
    (ROOT/'results/well_65_70.svg').write_text('\n'.join(p)+'\n',encoding='utf-8',newline='\n')
if __name__=='__main__':
    if sys.argv[1]=='freeze':freeze()
    elif sys.argv[1]=='evaluate':evaluate()
    else:raise ValueError('Choose freeze or evaluate')
