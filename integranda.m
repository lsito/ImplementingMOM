function f = integranda(delta_z, k, a, z_m, z, x)
    f = (4*pi)*exp(-1i*k*sqrt(a^2 + (z_m - (x*delta_z + z))^2))/sqrt(a^2 + (z_m - (x*delta_z + z))^2);
end
