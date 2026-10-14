"""Reproduce population and preregistered Fourier analyses using Python stdlib.

Run `python analyze.py` from any directory. No Gaussian or packages required.
Input CSVs contain audited extracted values, not synthetic energies.
"""
from pathlib import Path
import csv,math,json
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'
R=.00831446261815324
T=298.
HA=2625.4996394799
def read(name):
    return list(csv.DictReader((ROOT/'data'/name).open(encoding='utf-8-sig')))
def write(name,rows):
    with (OUT/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def solve(a,b):
    m=[list(row)+[v] for row,v in zip(a,b)];n=len(b)
    for k in range(n):
        j=max(range(k,n),key=lambda j:abs(m[j][k]));m[k],m[j]=m[j],m[k]
        if abs(m[k][k])<1e-12:raise ValueError('Singular model matrix')
        v=m[k][k];m[k]=[x/v for x in m[k]]
        for j in range(n):
            if j!=k:
                v=m[j][k];m[j]=[x-v*y for x,y in zip(m[j],m[k])]
    return [row[-1] for row in m]
def features(deg,n):return [math.cos(math.radians(deg)*k) for k in range(n+1)]
def value(deg,c):return sum(x*y for x,y in zip(features(deg,len(c)-1),c))
def fit(rows,n):
    x=[features(float(r['target_angle_deg']),n) for r in rows]
    y=[float(r['relative_electronic_kj_mol']) for r in rows]
    return solve([[sum(row[i]*row[j] for row in x) for j in range(n+1)] for i in range(n+1)], [sum(row[i]*v for row,v in zip(x,y)) for i in range(n+1)])
def roots(c):
    def derivative(a):return -sum(k*c[k]*math.sin(math.radians(a)*k) for k in range(1,len(c)))
    found=[]
    for i in range(1800):
        a=i/10;b=(i+1)/10
        if derivative(a)*derivative(b)<0:
            for _ in range(45):
                mid=(a+b)/2
                if derivative(a)*derivative(mid)<=0:b=mid
                else:a=mid
            q=(a+b)/2
            curvature=-sum(k*k*c[k]*math.cos(math.radians(q)*k) for k in range(1,len(c)))
            found.append((q,'minimum' if curvature>0 else 'maximum'))
    return found
def integrate(c,boundary,n):
    # Midpoint quadrature on 0..180; mirror symmetry supplies the other half.
    # Integrate basins separately so the boundary never cuts a quadrature cell.
    def integ(lo,hi):return (hi-lo)/n*sum(math.exp(-value(lo+(i+.5)*(hi-lo)/n,c)/(R*T)) for i in range(n))
    gauche=integ(0,boundary);anti=integ(boundary,180)
    return 100*gauche/(gauche+anti)
def svg(rows,models,filename='fourier.svg'):
    def X(x):return 65+x*3.6
    def Y(y):return 390-y*12
    parts=['<svg xmlns="http://www.w3.org/2000/svg" width="800" height="500" viewBox="0 0 800 500">','<rect width="800" height="500" fill="white"/>','<g font-family="sans-serif" font-size="14" fill="#182536">','<text x="65" y="28" font-size="20">Butane: calculated points and Fourier models</text>']
    for y in range(0,26,5):parts += [f'<path d="M65 {Y(y)} H715" stroke="#dde3e9"/>',f'<text x="30" y="{Y(y)+5}">{y}</text>']
    for x in range(0,181,30):parts.append(f'<text x="{X(x)-8}" y="415">{x}</text>')
    for (n,c),color in zip(models.items(),['#2864b4','#bd6114']):
        points=' '.join(f'{X(a):.2f},{Y(value(a,c)):.2f}' for a in range(181))
        parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2" stroke-dasharray="6 4"/>')
        parts.append(f'<text x="{65+(3-n)*330}" y="475" fill="{color}">Dashed: cosine order {n} prediction</text>')
    for r in rows:
        x=X(float(r['target_angle_deg']));y=Y(float(r['relative_electronic_kj_mol']))
        if r['role']=='train':parts.append(f'<circle cx="{x}" cy="{y}" r="5" fill="#182536"/>')
        else:parts.append(f'<rect x="{x-5}" y="{y-5}" width="10" height="10" fill="#ffffff" stroke="#a32954" stroke-width="2"/>')
    parts+=['<text x="65" y="445">Angle (degrees); energy relative to 180 degrees (kJ/mol)</text>','<text x="65" y="55">Filled circles: 4 fit points | Open squares: 3 held-out calculations</text>','</g></svg>']
    (OUT/filename).write_text('\n'.join(parts)+'\n',encoding='utf-8')
def main():
    OUT.mkdir(exist_ok=True)
    rows=read('torsion_audit.csv');train=[r for r in rows if r['role']=='train'];test=[r for r in rows if r['role']=='heldout']
    models={n:fit(train,n) for n in (3,2)}
    residuals=[];metrics=[];basins=[]
    for n,c in models.items():
        errors=[]
        for r in test:
            a=float(r['target_angle_deg']);truth=float(r['relative_electronic_kj_mol']);pred=value(a,c);err=pred-truth;errors.append(err)
            residuals.append(dict(order=n,angle_deg=a,calculated_kj_mol=truth,predicted_kj_mol=pred,error_kj_mol=err))
        metrics.append(dict(order=n,n_train=len(train),n_heldout=len(test),mae_kj_mol=sum(map(abs,errors))/len(errors),rmse_kj_mol=math.sqrt(sum(e*e for e in errors)/len(errors)),max_abs_kj_mol=max(map(abs,errors))))
        stationary=roots(c);boundaries=[a for a,kind in stationary if kind=='maximum']
        if len(boundaries)==1:
            for grid in (1000,10000):basins.append(dict(order=n,boundary_deg=boundaries[0],points_per_basin=grid,total_gauche_percent=integrate(c,boundaries[0],grid),measure='uniform dphi; mirrored basins; electronic model only'))
        else:basins.append(dict(order=n,boundary_deg='',points_per_basin='',total_gauche_percent='',measure='No unique interior maximum; no three-basin interpretation'))
    write('heldout_residuals.csv',residuals);write('validation_metrics.csv',metrics);write('one_dimensional_populations.csv',basins)
    (OUT/'coefficients.json').write_text(json.dumps({str(n):{'coefficients_kj_mol':c,'stationary_points_0_180':roots(c)} for n,c in models.items()},indent=2)+'\n',encoding='utf-8')
    minima=read('minima_audit.csv');ref=next(r for r in minima if r['state']=='anti');g=next(r for r in minima if r['state']=='gauche_plus')
    pop=[]
    for name,gap,count in [('reported_gap_one_gauche',4.07,1),('reported_gap_two_gauche',4.07,2),('minimum_electronic', (float(g['optimization_energy_hartree'])-float(ref['optimization_energy_hartree']))*HA,2),('minimum_RRHO_G',(float(g['electronic_plus_gibbs_hartree'])-float(ref['electronic_plus_gibbs_hartree']))*HA,2)]:
        weight=math.exp(-gap/(R*T));total=count*weight/(1+count*weight)
        pop.append(dict(model=name,gap_kj_mol=gap,temperature_K=T,gauche_degeneracy=count,anti_percent=100*(1-total),one_gauche_percent=100*total/count,total_gauche_percent=100*total))
    write('populations.csv',pop);svg(rows,models);svg(rows,{},'calculated_points.svg')
    print(json.dumps({'metrics':metrics,'basin_models':basins,'population_models':pop},indent=2))
if __name__=='__main__':main()
