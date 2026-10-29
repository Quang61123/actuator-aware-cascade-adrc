function f = mav_frame(msgid, payload, crcExtra, seq)
%MAV_FRAME MAVLink 2 frame (no signing) with trailing-zero payload truncation.
%   System id 1, component id 51 (simulator).
persistent tab
if isempty(tab), tab = mav_crc_table(); end
pl = uint8(payload(:)');
n = numel(pl); while n > 1 && pl(n) == 0, n = n - 1; end
pl = pl(1:n);
hdr = uint8([253, n, 0, 0, mod(seq,256), 1, 51, bitand(msgid,255), ...
             bitand(bitshift(msgid,-8),255), bitshift(msgid,-16)]);
crc = mav_crc([hdr(2:end), pl, uint8(crcExtra)], tab);
f = [hdr, pl, uint8([bitand(crc,255), bitshift(crc,-8)])];
end

function crc = mav_crc(b, tab)
crc = 65535;
for i = 1:numel(b)
    crc = bitxor(bitshift(crc,-8), tab(bitxor(bitand(crc,255), double(b(i))) + 1));
end
end

function tab = mav_crc_table()
% table for CRC-16/MCRF4XX (X.25, reflected polynomial 0x8408)
tab = zeros(1,256);
for v = 0:255
    c = v;
    for j = 1:8
        if bitand(c,1), c = bitxor(bitshift(c,-1), 33800); else, c = bitshift(c,-1); end
    end
    tab(v+1) = c;
end
end
