function s = prof_pulse(t, period, phase)
    s = 1 * (mod(t-phase*0.5*period/pi, period)/period < 0.5);
end