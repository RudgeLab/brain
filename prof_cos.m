function s = prof_cos(t, period, phase)
    s = 0.5 * (1  + cos(2 * pi * t / period + phase));
end