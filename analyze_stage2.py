"""Freeze stage-2 models before jobs, then evaluate new reference values."""
import sys,json,math,hashlib
from analyze import ROOT,read,write,fit,value,roots,integrate
FROZEN=ROOT/'data/stage2_frozen_models.json'
ANGLES=[20,75,105,135,165]
def development():
    return read('torsion_audit.csv')+read('followup_15_45_audit.csv')
def freeze():
    rows=development()
    record={'date':'2026-11-11','fit_angles':sorted(float(r['target_angle_deg']) for r in rows),'validation_angles':ANGLES,
      'source_hashes':{n:hashlib.sha256((ROOT/'data'/n).read_bytes()).hexdigest() for n in ['torsion_audit.csv','followup_15_45_audit.csv']},
      'models':{str(n):{'role':role,'coefficients_kj_mol':fit(rows,n)} for n,role in [(3,'baseline'),(4,'primary'),(5,'secondary')]}}
    for m in record['models'].values():m['prospective_predictions']={str(a):value(a,m['coefficients_kj_mol']) for a in ANGLES}
    text=json.dumps(record,indent=2)+'\n'
    if FROZEN.exists() and FROZEN.read_text()!=text:raise ValueError('Refuse to change frozen models')
    FROZEN.write_text(text,encoding='utf-8',newline='\n')
    print('Frozen coefficients and predictions:',FROZEN.name)
def evaluate():
    record=json.loads(FROZEN.read_text());rows=read('stage2_audit.csv');res=[];metrics=[];diagnostics=[]
    assert sorted(int(r['target_angle_deg']) for r in rows)==ANGLES
    for name,sha in record['source_hashes'].items():
        assert hashlib.sha256((ROOT/'data'/name).read_bytes()).hexdigest()==sha
    for order,m in record['models'].items():
        c=m['coefficients_kj_mol'];errors=[]
        for r in rows:
            a=int(r['target_angle_deg']);truth=float(r['relative_electronic_kj_mol']);pred=m['prospective_predictions'][str(a)];err=pred-truth;errors.append(err)
            res.append(dict(order=order,role=m['role'],angle_deg=a,reference_kj_mol=truth,prediction_kj_mol=pred,error_kj_mol=err))
        metrics.append(dict(order=order,role=m['role'],n_fit=9,n_validation=5,mae_kj_mol=sum(map(abs,errors))/5,rmse_kj_mol=math.sqrt(sum(e*e for e in errors)/5),max_abs_kj_mol=max(map(abs,errors))))
        stationary=roots(c)
        good=len(stationary)==2 and [kind for _,kind in stationary]==['minimum','maximum']
        # Check endpoint curvatures too: 0 is the syn barrier; 180 the anti minimum.
        good=good and -sum(k*k*c[k] for k in range(1,len(c)))<0 and -sum(k*k*c[k]*(-1)**k for k in range(1,len(c)))>0
        for grid in [1000,10000]:
            diagnostics.append(dict(order=order,expected_topology=bool(good),gauche_minimum_deg=stationary[0][0] if good else '',gauche_gap_kj_mol=value(stationary[0][0],c)-value(180,c) if good else '',boundary_deg=stationary[1][0] if good else '',cells_per_basin=grid,total_gauche_percent=integrate(c,stationary[1][0],grid) if good else ''))
    write('stage2_residuals.csv',res);write('stage2_metrics.csv',metrics);write('stage2_basin_diagnostics.csv',diagnostics)
    plot(record,rows)
    print(json.dumps({'metrics':metrics,'diagnostics':diagnostics},indent=2))
def plot(record,validation):
    def x(a):return 70+3.5*a
    def y(e):return 380-12*e
    p=['<svg xmlns="http://www.w3.org/2000/svg" width="820" height="510" viewBox="0 0 820 510">','<rect width="820" height="510" fill="white"/>','<g font-family="sans-serif" font-size="14">','<text x="70" y="25" font-size="21">Stage 2: nine fit points and five new calculations</text>']
    for e in range(0,26,5):p += [f'<path d="M70 {y(e)} H700" stroke="#ddd"/>',f'<text x="35" y="{y(e)+5}">{e}</text>']
    for a in range(0,181,30):p.append(f'<text x="{x(a)-8}" y="405">{a}</text>')
    for i,(order,m) in enumerate(record['models'].items()):
        color=['#777777','#176cc0','#c66a18'][i];c=m['coefficients_kj_mol']
        pts=' '.join(f'{x(a):.2f},{y(value(a,c)):.2f}' for a in range(181))
        p += [f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2" stroke-dasharray="5 4"/>',f'<text x="{70+i*230}" y="475" fill="{color}">Order {order}: {m["role"]}</text>']
    for r in development():p.append(f'<circle cx="{x(float(r["target_angle_deg"]))}" cy="{y(float(r["relative_electronic_kj_mol"]))}" r="4" fill="#222"/>')
    for r in validation:p.append(f'<rect x="{x(float(r["target_angle_deg"]))-5}" y="{y(float(r["relative_electronic_kj_mol"]))-5}" width="10" height="10" fill="white" stroke="#ac204a" stroke-width="2"/>')
    p+=['<text x="70" y="50">Filled: development | Open squares: new validation | Dashed: models</text>','<text x="70" y="440">Angle (degrees); relative electronic energy (kJ/mol)</text>','</g></svg>']
    (ROOT/'results/stage2.svg').write_text('\n'.join(p)+'\n',encoding='utf-8',newline='\n')
if __name__=='__main__':
    if sys.argv[1]=='freeze':freeze()
    elif sys.argv[1]=='evaluate':evaluate()
    else:raise ValueError('Choose freeze or evaluate')
