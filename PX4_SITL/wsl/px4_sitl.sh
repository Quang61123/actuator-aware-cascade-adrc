#!/bin/bash
# PX4 SITL control script with instance support (WSL distribution PX4-Ubuntu-22.04).
#   px4_sitl.sh start <ctrlId> <scen> <tag> [inst] [src]   relay + PX4 instance inst (default 0);
#                                                    src 1 = controller reads the raw measurement
#   px4_sitl.sh arm [inst]                           commander arm -f
#   px4_sitl.sh stop [inst]                          kill this instance and its relay
#   px4_sitl.sh lastlog [inst]                       newest ulog of the instance
# Instance i: PX4 simulator TCP port 4560+i, MATLAB relay port 4600+i,
# airframe 1011i (copy of 10099_ladrc_sitl with the run's controller/scenario),
# working directory build/px4_sitl_default/instance_i.
B=$HOME/myPX4/PX4-Autopilot/build/px4_sitl_default
W=${PX4SITL_WORK_WSL:-$(cd "$(dirname "$0")" && pwd)}   # work folder (where this script is copied by sitl_run)
case "$1" in
start)
  I=${5:-0}; AFID=$((10110+I)); WD=$B/instance_$I
  pkill -9 -f "bin/px4 -i $I " 2>/dev/null; pkill -9 -f "sim_relay.py $((4600+I)) " 2>/dev/null
  for k in 1 2 3 4 5 6 7 8 9 10; do pgrep -f "bin/px4 -i $I " >/dev/null || break; sleep 0.3; done
  rm -f /tmp/px4_lock-$I /tmp/px4-sock-$I
  rm -rf $WD; mkdir -p $WD
  rm -f $B/etc/init.d-posix/airframes/${AFID}_*
  sed -e "s/@CTRL@/$2/" -e "s/@SCEN@/$3/" -e "s/@SRC@/${6:-0}/" $W/10099_ladrc_sitl > $B/etc/init.d-posix/airframes/${AFID}_ladrc_sitl_i$I
  setsid nohup python3 $W/sim_relay.py $((4600+I)) $((4560+I)) > $W/relay_$4.log 2>&1 < /dev/null &
  sleep 0.5
  cd $WD
  PX4_SYS_AUTOSTART=$AFID setsid nohup $B/bin/px4 -i $I -d $B/etc > $W/px4_$4.log 2>&1 < /dev/null &
  sleep 0.5
  echo started ;;
arm)
  I=${2:-0}; cd $B/instance_$I && $B/bin/px4-commander --instance $I arm -f ;;
stop)
  I=${2:-0}
  pkill -9 -f "bin/px4 -i $I " 2>/dev/null; pkill -9 -f "sim_relay.py $((4600+I)) " 2>/dev/null; sleep 0.3
  rm -f /tmp/px4_lock-$I /tmp/px4-sock-$I; echo stopped ;;
lastlog)
  I=${2:-0}; ls -t $(find $B/instance_$I/log -name "*.ulg" 2>/dev/null) 2>/dev/null | head -1 ;;
esac
