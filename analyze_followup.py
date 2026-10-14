"""Evaluate 15/45 degree targeted checks using frozen original coefficients."""
import json,math
from analyze import ROOT,read,write,value,svg
def main():
    frozen=json.loads((ROOT/'data/frozen_models_3090f6b.json').read_text())
    rows=read('followup_15_45_audit.csv');residuals=[];metrics=[]
    for order,model in frozen.items():
        c=model['coefficients_kj_mol'];errors=[]
        for r in rows:
            a=float(r['target_angle_deg']);truth=float(r['relative_electronic_kj_mol']);pred=value(a,c);err=pred-truth;errors.append(err)
            residuals.append(dict(order=order,angle_deg=a,calculated_kj_mol=truth,predicted_kj_mol=pred,error_kj_mol=err,selection='targeted after inspecting original 30-degree error; no refit'))
        metrics.append(dict(order=order,n_targeted_followup=len(rows),mae_kj_mol=sum(map(abs,errors))/len(errors),rmse_kj_mol=math.sqrt(sum(e*e for e in errors)/len(errors)),max_abs_kj_mol=max(map(abs,errors))))
    write('followup_residuals.csv',residuals);write('followup_metrics.csv',metrics)
    # Original seven-point figures remain unchanged. This plot labels both
    # original held-out and targeted follow-up points as evaluation points.
    svg(read('torsion_audit.csv')+rows,{int(n):m['coefficients_kj_mol'] for n,m in frozen.items()},'followup_nine_points.svg')
    p=ROOT/'results/followup_nine_points.svg'
    p.write_text(p.read_text().replace('3 held-out calculations','5 evaluation calculations (3 original + 2 targeted)'),encoding='utf-8',newline='\n')
    print(json.dumps({'residuals':residuals,'metrics':metrics},indent=2))
if __name__=='__main__':main()
