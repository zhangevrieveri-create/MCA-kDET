"""Create an estimated-only confirmation entry; original pilot sources stay intact."""
from pathlib import Path
ROOT=Path(__file__).resolve().parent


def once(text,old,new):
    assert text.count(old)==1,(old,text.count(old))
    return text.replace(old,new)


def main():
    s=(ROOT/'sampling_calibration.R').read_text(encoding='utf-8')
    s=once(s,'for(phase in c("oracle","estimated"))','for(phase in "estimated")')
    s=once(s,'paths <- setNames(file.path(config$outdir,c("oracle","estimated"),"chunks",filename),c("oracle","estimated"))',
             'paths <- setNames(file.path(config$outdir,"estimated","chunks",filename),"estimated")')
    s=once(s,'rows <- list(oracle=list(),estimated=list())','rows <- list(estimated=list())')
    s=once(s,'keep=r==1L','keep=(r==1L || r==101L)')
    s=once(s,'if(r==1L) saveRDS','if(r==1L || r==101L) saveRDS')
    s=once(s,'paste0(job$scenario,"_r00001.rds")','sprintf("%s_r%05d.rds",job$scenario,r)')
    s=once(s,'scenario=scenario,rep=r,d=config$d','scenario=scenario,rep=r,seed=sampling_seed(case,r,0L),d=config$d')
    (ROOT/'sampling_calibration_confirmation.R').write_text(s,encoding='utf-8')
    d=(ROOT/'run_sampling_validation_fullfit.R').read_text(encoding='utf-8')
    d=d.replace('sampling_calibration.R','sampling_calibration_confirmation.R')
    d=d.replace('run_sampling_validation_fullfit.R','run_sampling_confirmation.R')
    d=d.replace('sampling_fullfit_pilot_d5','sampling_fullfit_confirmation_d5')
    d=once(d,'config <- list(d=5L','config <- list(phase="estimated",d=5L')
    d=once(d,'arg("reps","100")','arg("reps","1000")')
    d=once(d,'arg("first","1")','arg("first","101")')
    d=d.replace('c("oracle","estimated")','"estimated"')
    (ROOT/'run_sampling_confirmation.R').write_text(d,encoding='utf-8')
    print('Created estimated-only confirmation entry; source pilot files unchanged.')


if __name__=='__main__':main()
