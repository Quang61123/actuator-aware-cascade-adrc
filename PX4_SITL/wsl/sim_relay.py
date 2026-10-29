#!/usr/bin/env python3
"""TCP relay inside WSL between PX4 SITL (simulator_mavlink, client on port 4560)
and the MATLAB plant on Windows (client on port 4561, reached through WSL
localhost forwarding). Needed because the Windows firewall blocks inbound
connections from WSL to MATLAB and the account cannot add a firewall rule.
Bytes are forwarded unchanged in both directions (TCP_NODELAY)."""
import select, socket, sys

def listen(host, port):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    s.bind((host, port)); s.listen(1)
    return s

def main():
    pm = int(sys.argv[1]) if len(sys.argv) > 1 else 4561
    pp = int(sys.argv[2]) if len(sys.argv) > 2 else 4560
    lm = listen('0.0.0.0', pm)
    print('relay: waiting for MATLAB on', pm, flush=True)
    m, _ = lm.accept(); m.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
    lp = listen('127.0.0.1', pp)
    print('relay: MATLAB connected; waiting for PX4 on', pp, flush=True)
    p, _ = lp.accept(); p.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
    print('relay: PX4 connected', flush=True)
    peer = {m: p, p: m}
    while True:
        r, _, _ = select.select([m, p], [], [])
        for s in r:
            d = s.recv(65536)
            if not d:
                print('relay: connection closed', flush=True); return
            peer[s].sendall(d)

if __name__ == '__main__':
    main()
