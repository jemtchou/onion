mkdir -p jobs
cd jobs
rm -f *
begin=$1
end=$(($2-1))
for f in `seq $1 $end`
do
  echo "Job $f"
  cat << 'EOF' > run_$f.sh
  #PBS -N geant4_%GEO%_%RUN% 
  #PBS -l walltime=23:50:00
  #PBS -o lxpub01.jinr.ru:/afs/jinr.ru/user/j/jemtchou/onion/prod
  cd $TMPDIR
  cp /afs/jinr.ru/user/j/jemtchou/onion/run/run_template.mac ./init.mac
  source /afs/.jinr.ru/geant4/setup_geant4.10.06.sh
  /afs/jinr.ru/user/j/jemtchou/onion/build/onion init.mac %RUN% >& run_%RUN%.out
  cat ./init.mac >> run_%RUN%.log
#  mv pos.csv pos_%RUN%.csv
  ls -l
  cp fout_0.csv /afs/jinr.ru/user/j/jemtchou/onion/prod/%RUN%.csv
  cp run_%RUN%.out /afs/jinr.ru/user/j/jemtchou/onion/prod/
EOF
  sed -i "s/%RUN%/$f/g" run_$f.sh
  sed -i "s/%GEO%/$geo/g" run_$f.sh
 
  qsub -q common run_$f.sh
done

