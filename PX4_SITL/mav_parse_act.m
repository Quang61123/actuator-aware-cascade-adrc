function [act, buf] = mav_parse_act(buf)
%MAV_PARSE_ACT Extract HIL_ACTUATOR_CONTROLS (id 93) messages from a byte
%   buffer (MAVLink 1 or 2). Other messages are skipped. Returns a struct
%   array with fields t (time_usec), u (16 controls), mode; and the
%   unconsumed remainder of the buffer. CRC is not re-checked (TCP stream).
act = struct('t', {}, 'u', {}, 'mode', {});
i = 1; n = numel(buf);
while i <= n
    b = buf(i);
    if b == 253
        if n - i + 1 < 10, break; end
        len = double(buf(i+1)); inc = buf(i+2);
        tot = 12 + len + 13*double(bitand(inc,1));
        if n - i + 1 < tot, break; end
        id = double(buf(i+7)) + 256*double(buf(i+8)) + 65536*double(buf(i+9));
        if id == 93
            pl = zeros(1,81,'uint8'); m = min(len,81); pl(1:m) = buf(i+10:i+9+m);
            act(end+1) = struct('t', double(typecast(pl(1:8),'uint64')), ...
                'u', double(typecast(pl(17:80),'single')), 'mode', double(pl(81))); %#ok<AGROW>
        end
        i = i + tot;
    elseif b == 254
        if n - i + 1 < 6, break; end
        len = double(buf(i+1)); tot = 8 + len;
        if n - i + 1 < tot, break; end
        if buf(i+5) == 93
            pl = zeros(1,81,'uint8'); m = min(len,81); pl(1:m) = buf(i+6:i+5+m);
            act(end+1) = struct('t', double(typecast(pl(1:8),'uint64')), ...
                'u', double(typecast(pl(17:80),'single')), 'mode', double(pl(81))); %#ok<AGROW>
        end
        i = i + tot;
    else
        i = i + 1;
    end
end
buf = buf(i:end);
end
