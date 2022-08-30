#USAGE= gawk -f pro2log_ini.awk Background_P_file(no framed), Observed_P_file(framed)
BEGIN{
  while(getline<ARGV[1]>0) # read backgrounb probabilities
    BP[$1,$2]=$3;

  ARGV[1]="";
}
{
  print $1, $2, $3, log($4/BP[$2,$3]);
}
