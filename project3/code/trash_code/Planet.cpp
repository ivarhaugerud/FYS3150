

Planet(m, x0, y0, z0, vx0, vy0, vz0);
{
  mass = m;

  initial_position(0) =x0;
  initial_position(1) =y0;
  initial_position(2) =z0;

  intitial_velocity(0) = vx0;
  intitial_velocity(1) = vy0;
  intitial_velocity(2) = vz0;

}

void change_initial_position(x_new, y_new, z_new)
{
  initial_position(0) = x_new;
  initial_position(1) = y_new;
  initial_position(2) = z_new;
}

void change_initial_velocity(vx_new, vy_new, vz_new)
{
  intitial_velocity(0) = vx_new;
  intitial_velocity(1) = vy_new;
  intitial_velocity(2) = vz_new;
}

void change_mass(m_new)
{
  mass = m_new;
}
