#!/bin/bash
# Audit 297: czas do zamkniecia kanalu E1, para/orto x zachowawcza/retardowana,
# ziarna 40..49 sparowane.  Po DWA procesy naraz (trzeci rdzen zajmuje dlugi
# przebieg z 296g, czwarty zostaje wolny -- budzet jest zegarowy).
S="$1"; D="$S/close"
jobs_list=()
for cfg in zach ret; do for ph in 1 2; do for seed in $(seq 40 49); do
  jobs_list+=("$cfg $ph $seed"); done; done; done
run_one() {
  set -- $1
  cfg=$1; ph=$2; seed=$3
  envs="CREM_ACTION_TRACE=1"
  [ "$cfg" = ret ] && envs="$envs CREM_RETARDED_BETWEEN_PHOTONS=1 CREM_PRESCRIBED_EMISSION_PATTERN=1"
  env $envs CREM_PROBE_SEED=$seed timeout 260 "$S/close/stopcause" 1 220 $ph \
     2> "$D/$cfg-$ph-$seed.err" > "$D/$cfg-$ph-$seed.out"
}
i=0
for j in "${jobs_list[@]}"; do
  run_one "$j" &
  i=$((i+1)); if [ $((i%2)) -eq 0 ]; then wait; fi
done
wait
echo DONE > "$D/DONE"
