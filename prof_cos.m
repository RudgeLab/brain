function s = prof_cos(t, period, phase) % function to generate a cosine wave
    s = 0.5 * (1  + cos(2 * pi * t / period + phase));
end