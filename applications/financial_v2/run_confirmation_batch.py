"""Run the immediate confirmation job, log progress, and resume R checkpoints."""
from pathlib import Path
import argparse
import json
import os
import re
import subprocess
import sys
import threading
import time
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parent
R=Path('D:/R-4.6.1/bin/Rscript.exe')
CASES=['null_iid','null_garch95','null_garch99','null_heterogeneous','null_garch95_t5',
       'null_garch999','alt_sector_shape','alt_one_bank','null_equal_mixtures']


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--reps',type=int,default=1000);p.add_argument('--first',type=int,default=101)
    p.add_argument('--B',type=int,default=999);p.add_argument('--cores',type=int,default=30)
    p.add_argument('--out',type=Path,default=ROOT/'sampling_fullfit_confirmation_d5')
    p.add_argument('--dry-run',action='store_true');a=p.parse_args()
    out=a.out.resolve()
    command=[str(R),'--vanilla',(ROOT/'run_sampling_confirmation.R').as_posix(),
             f'--reps={a.reps}',f'--first={a.first}',f'--B={a.B}',f'--cores={a.cores}',
             '--scenarios='+','.join(CASES),'--out='+out.as_posix()]
    if a.dry_run:
        print(json.dumps(dict(command=command,primary_stride=4,scientific_status='candidate awaiting independent confirmation'),ensure_ascii=False));return
    out.mkdir(parents=True,exist_ok=True)
    logs=out/'execution_logs';logs.mkdir(exist_ok=True)
    env=os.environ.copy()
    for key in ('LC_ALL','LC_CTYPE','LANG'):env.pop(key,None)
    env['R_LIBS_USER']='C:/Users/zhang/AppData/Local/R/win-library/4.6'
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):env[key]='1'
    (out/'execution_request.json').write_text(json.dumps(dict(command=command,primary_stride=4,
        first_rep=a.first,last_rep=a.first+a.reps-1,replications_per_case=a.reps,B=a.B,cores=a.cores,
        source='Independent confirmation uses new simulated panels; real p-values not computed'),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    started=time.monotonic();last_print=None
    def progress(stage,pid=None):
        nonlocal last_print
        done={case:set() for case in CASES}
        for file in (out/'estimated/chunks').glob('*.csv'):
            match=re.fullmatch(r'(.+)__r(\d+)-(\d+)\.csv',file.name)
            if match and match[1] in done:
                done[match[1]].update(r for r in range(int(match[2]),int(match[3])+1) if a.first<=r<a.first+a.reps)
        counts={case:len(values) for case,values in done.items()}
        result=dict(updated_utc=datetime.now(timezone.utc).isoformat(),stage=stage,python_pid=os.getpid(),R_launcher_pid=pid,
                    estimated_replications_completed=counts,total_completed=sum(counts.values()),total_planned=a.reps*len(CASES),
                    elapsed_minutes=round((time.monotonic()-started)/60,2),primary_stride=4,
                    scientific_status='Not yet reviewed; numerical completion does not establish validity',real_pvalues_computed=False)
        (out/'batch_progress.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        signature=(stage,result['total_completed'])
        if signature!=last_print:
            print(f"[{stage}] {result['total_completed']}/{result['total_planned']} independent outer repetitions; elapsed {result['elapsed_minutes']:.1f} min",flush=True)
            last_print=signature
        return result
    for attempt in range(1,4):
        logfile=logs/f'R_attempt_{attempt:02d}.log'
        with logfile.open('w',encoding='utf-8') as log:
            proc=subprocess.Popen(command,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                  text=True,encoding='utf-8',errors='replace',bufsize=1)
            def consume():
                for line in proc.stdout:
                    log.write(line);log.flush();print(line.rstrip(),flush=True)
            thread=threading.Thread(target=consume,daemon=True);thread.start()
            progress('RUNNING',proc.pid)
            while True:
                try:code=proc.wait(timeout=30);break
                except subprocess.TimeoutExpired:progress('RUNNING',proc.pid)
            thread.join(timeout=10)
        if code==0:break
        progress(f'R_EXIT_{code}_RETRY_{attempt}',proc.pid)
    else:
        progress('FAILED_NEEDS_REVIEW');raise SystemExit('R failed three times; all completed chunks retained. Inspect execution_logs.')
    status=progress('SUMMARIZING')
    if any(n!=a.reps for n in status['estimated_replications_completed'].values()):
        progress('INCOMPLETE_NEEDS_REVIEW');raise SystemExit('R exited but expected confirmation replications are incomplete')
    summary=subprocess.run([sys.executable,str(ROOT/'summarize_validation.py'),str(out/'estimated'),'--main-length=4'],
                           env=env,capture_output=True,text=True,encoding='utf-8',errors='replace')
    (logs/'summary.log').write_text(summary.stdout+summary.stderr,encoding='utf-8')
    if summary.returncode:progress('SUMMARY_FAILED_NEEDS_REVIEW');raise SystemExit(summary.stderr)
    progress('COMPLETED_AWAITING_SCIENTIFIC_REVIEW')
    print('Confirmation complete; inspect estimated/calibration_rates.csv, calibration_fwer.csv and the Chinese report. No real-data p-values or paper edits.',flush=True)


if __name__=='__main__':main()
