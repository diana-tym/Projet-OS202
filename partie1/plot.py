import matplotlib.pyplot as plt

# Données : nombre de threads et temps mesurés pour chaque métrique
threads = [1, 2, 4, 6, 8, 16]

# Temps total de la simulation en millisecondes (baseline pour 1 thread = 959 ms)
sim_total = [959, 943, 894, 837, 953, 1823]

# Temps moyen pour la mise à jour du modèle en microsecondes (baseline pour 1 thread = 103.316 µs)
update_times = [103.316, 191.909, 195.352, 361.533, 680.83, 1889.6]

# Temps moyen pour l'affichage en microsecondes (baseline pour 1 thread = 1574.56 µs)
display_times = [1574.56, 1460.3, 1369.05, 1107.66, 990.679, 1327.66]

# Calcul des speedups : speedup = temps_baseline / temps_actuel
sim_speedup = [sim_total[0] / t for t in sim_total]
update_speedup = [update_times[0] / t for t in update_times]
display_speedup = [display_times[0] / t for t in display_times]

# Affichage des courbes
plt.figure(figsize=(10, 6))
plt.plot(threads, sim_speedup, marker='o', label='Speedup Simulation Total')
plt.plot(threads, update_speedup, marker='s', label='Speedup Model Update')
plt.plot(threads, display_speedup, marker='^', label='Speedup Display')
plt.xlabel('Nombre de Threads')
plt.ylabel('Speedup (Baseline: 1 thread)')
plt.title('Speedup en fonction du nombre de Threads')
plt.legend()
plt.grid(True)
plt.xticks(threads)  # Pour afficher les valeurs exactes des threads
plt.show()
