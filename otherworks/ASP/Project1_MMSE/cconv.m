function val=cconv(a,b)

pad=length(a)-length(b);

if pad > 0,
    val=ifft(fft(a).*fft([b zeros(1,abs(pad))]));
elseif pad < 0,
    val=ifft(fft(b).*fft([a zeros(1,abs(pad))]));
else
    val=ifft(fft(a).*fft(b));
end