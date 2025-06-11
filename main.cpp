#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <limits>
#include <iomanip>


using ProcessingTimes = std::vector<std::vector<double>>;
using Permutation = std::vector<int>;

// Funkcja do obliczania Cmax dla danej permutacji
// processing_times[job_idx][machine_idx]
// permutation[i] = job_idx
double calculate_cmax(const ProcessingTimes& p_times, const Permutation& perm) {
    if (perm.empty()) {
        return 0.0;
    }
    if (p_times.empty()) {
        return 0.0;
    }

    int n_jobs_in_perm = perm.size();
    int m_machines = p_times[0].size();

    if (m_machines == 0) {
        return 0.0;
    }

    // completion_times[i][j] = czas zakończenia i-tego zadania w permutacji na maszynie j
    std::vector<std::vector<double>> C(n_jobs_in_perm, std::vector<double>(m_machines, 0.0));

    for (int i = 0; i < n_jobs_in_perm; ++i) { // Iteracja po zadaniach w permutacji
        int current_job_idx = perm[i]; // Aktualny indeks zadania
        for (int j = 0; j < m_machines; ++j) { // Iteracja po maszynach
            double prev_job_same_machine_completion = (i == 0) ? 0.0 : C[i - 1][j];
            double same_job_prev_machine_completion = (j == 0) ? 0.0 : C[i][j - 1];

            C[i][j] = std::max(prev_job_same_machine_completion, same_job_prev_machine_completion) + p_times[current_job_idx][j];
        }
    }
    return C[n_jobs_in_perm - 1][m_machines - 1];
}

// Funkcja pomocnicza do wyświetlania permutacji i Cmax
void print_solution(const std::string& algorithm_name, const Permutation& perm, double cmax) {
    std::cout << algorithm_name << ":\n";
    std::cout << "  Permutacja: ";
    for (size_t i = 0; i < perm.size(); ++i) {
        std::cout << perm[i] + 1 << (i == perm.size() - 1 ? "" : " -> "); // Wyświetlamy zadania 1-indeksowane
    }
    std::cout << "\n  Cmax: " << std::fixed << std::setprecision(2) << cmax << "\n\n";
}

// 1. Przegląd zupełny
Permutation brute_force_algorithm(const ProcessingTimes& p_times, double& best_cmax_val) {
    int n_jobs = p_times.size();
    if (n_jobs == 0) {
        best_cmax_val = 0.0;
        return {};
    }


    Permutation current_perm(n_jobs);
    std::iota(current_perm.begin(), current_perm.end(), 0); // Inicjalizacja: 0, 1, ..., n-1

    Permutation best_perm = current_perm;
    best_cmax_val = std::numeric_limits<double>::max();

    do {
        double cmax = calculate_cmax(p_times, current_perm);
        if (cmax < best_cmax_val) {
            best_cmax_val = cmax;
            best_perm = current_perm;
        }
    } while (std::next_permutation(current_perm.begin(), current_perm.end()));

    return best_perm;
}

// 2. Algorytm NEH
Permutation neh_algorithm(const ProcessingTimes& p_times, double& best_cmax_val) {
    int n_jobs = p_times.size();
    if (n_jobs == 0) {
        best_cmax_val = 0.0;
        return {};
    }
    int m_machines = p_times[0].size();
    if (m_machines == 0) {
        best_cmax_val = 0.0;
        return {};
    }

    // Krok 1: Oblicz sumę czasów przetwarzania dla każdego zadania
    std::vector<std::pair<double, int>> job_total_times(n_jobs);
    for (int i = 0; i < n_jobs; ++i) {
        double total_time = 0;
        for (int j = 0; j < m_machines; ++j) {
            total_time += p_times[i][j];
        }
        job_total_times[i] = {total_time, i}; // {suma czasów, indeks zadania}
    }

    // Krok 2: Posortuj zadania malejąco wg sumy czasów
    std::sort(job_total_times.rbegin(), job_total_times.rend());

    Permutation current_best_perm;

    // Krok 3 i 4: Iteracyjne dodawanie zadań
    for (int k = 0; k < n_jobs; ++k) {
        int job_to_insert = job_total_times[k].second;

        Permutation best_insertion_perm;
        double min_cmax_for_insertion = std::numeric_limits<double>::max();

        // Wypróbuj wstawienie zadania na każdą możliwą pozycję (k+1 pozycji)
        for (int pos = 0; pos <= k; ++pos) {
            Permutation temp_perm = current_best_perm;
            temp_perm.insert(temp_perm.begin() + pos, job_to_insert);

            double cmax = calculate_cmax(p_times, temp_perm);

            if (cmax < min_cmax_for_insertion) {
                min_cmax_for_insertion = cmax;
                best_insertion_perm = temp_perm;
            }
        }
        current_best_perm = best_insertion_perm;
        best_cmax_val = min_cmax_for_insertion;
    }


    return current_best_perm;
}


// Algorytm Johnsona dla m=2
Permutation johnson_algorithm(const ProcessingTimes& p_times, double& best_cmax_val) {
    int n_jobs = p_times.size();
    if (n_jobs == 0) {
        best_cmax_val = 0.0;
        return {};
    }


    struct JobInfo {
        int id;
        double p1, p2;
    };

    std::vector<JobInfo> jobs(n_jobs);
    for (int i = 0; i < n_jobs; ++i) {
        jobs[i] = {i, p_times[i][0], p_times[i][1]};
    }

    std::vector<int> group1, group2;
    for (const auto& job : jobs) {
        if (job.p1 < job.p2) {
            group1.push_back(job.id);
        } else {
            group2.push_back(job.id);
        }
    }

    // Sortuj grupę 1 rosnąco wg p1
    std::sort(group1.begin(), group1.end(), [&](int a, int b) {
        return p_times[a][0] < p_times[b][0];
    });

    // Sortuj grupę 2 malejąco wg p2
    std::sort(group2.begin(), group2.end(), [&](int a, int b) {
        return p_times[a][1] > p_times[b][1];
    });

    Permutation optimal_perm;
    optimal_perm.insert(optimal_perm.end(), group1.begin(), group1.end());
    optimal_perm.insert(optimal_perm.end(), group2.begin(), group2.end());

    best_cmax_val = calculate_cmax(p_times, optimal_perm);
    return optimal_perm;
}

// 3. FNEH (NEH z akceleracją)
// 4.0 FNEH – ulepszony NEH (Fast NEH)
// Działa jak NEH, ale po każdej iteracji wyjmuje ostatnio dodane zadanie,
// i sprawdza, czy po ponownym wstawieniu w inne miejsce Cmax będzie lepszy.
Permutation fneh_algorithm(const ProcessingTimes& p_times, double& best_cmax_val) {
    int n_jobs = p_times.size();
    int m_machines = p_times[0].size();
    if (n_jobs == 0 || m_machines == 0) {
        best_cmax_val = 0.0;
        return {};
    }

    // Krok 1: policz sumy czasów przetwarzania dla każdego zadania
    std::vector<std::pair<double, int>> job_total_times(n_jobs);
    for (int i = 0; i < n_jobs; ++i) {
        job_total_times[i].first = std::accumulate(p_times[i].begin(), p_times[i].end(), 0.0);
        job_total_times[i].second = i;
    }

    // Krok 2: sortuj malejąco po sumie czasów
    std::sort(job_total_times.rbegin(), job_total_times.rend());

    Permutation perm; // aktualna najlepsza permutacja
    best_cmax_val = 0.0;

    for (int k = 0; k < n_jobs; ++k) {
        int job_to_add = job_total_times[k].second;

        Permutation best_partial_perm;
        double best_partial_cmax = std::numeric_limits<double>::max();

        // Spróbuj dodać nowe zadanie na każdą pozycję
        for (int pos = 0; pos <= perm.size(); ++pos) {
            Permutation temp = perm;
            temp.insert(temp.begin() + pos, job_to_add);
            double cmax = calculate_cmax(p_times, temp);

            if (cmax < best_partial_cmax) {
                best_partial_cmax = cmax;
                best_partial_perm = temp;
            }
        }

        perm = best_partial_perm;
        best_cmax_val = best_partial_cmax;

        // Akceleracja FNEH: usuń ostatnio dodane i przetestuj ponownie
        int last_added_job = job_to_add;
        perm.erase(std::find(perm.begin(), perm.end(), last_added_job));

        Permutation fneh_best_perm;
        double fneh_best_cmax = std::numeric_limits<double>::max();

        for (int pos = 0; pos <= perm.size(); ++pos) {
            Permutation temp = perm;
            temp.insert(temp.begin() + pos, last_added_job);
            double cmax = calculate_cmax(p_times, temp);
            if (cmax < fneh_best_cmax) {
                fneh_best_cmax = cmax;
                fneh_best_perm = temp;
            }
        }

        perm = fneh_best_perm;
        best_cmax_val = fneh_best_cmax;
    }

    return perm;
}


// 4. Algorytm Podziału i Ograniczeń (Branch and Bound)

// Funkcja do obliczania dolnego ograniczenia (LB) dla węzła (częściowej permutacji)
// scheduled_perm: już uszeregowane zadania
// unscheduled_jobs_indices: indeksy zadań, które jeszcze nie są uszeregowane
double calculate_lower_bound(const ProcessingTimes& p_times,
                             const Permutation& scheduled_perm,
                             const std::vector<bool>& is_scheduled) {
    int n_total_jobs = p_times.size();
    if (n_total_jobs == 0) return 0.0;
    int m_machines = p_times[0].size();
    if (m_machines == 0) return 0.0;

    // Oblicz czasy zakończenia dla już uszeregowanych zadań
    std::vector<std::vector<double>> C_partial(scheduled_perm.size(), std::vector<double>(m_machines, 0.0));
    if (!scheduled_perm.empty()) {
        for (size_t i = 0; i < scheduled_perm.size(); ++i) {
            int current_job_idx = scheduled_perm[i];
            for (int j = 0; j < m_machines; ++j) {
                double prev_job_same_machine_completion = (i == 0) ? 0.0 : C_partial[i - 1][j];
                double same_job_prev_machine_completion = (j == 0) ? 0.0 : C_partial[i][j - 1];
                C_partial[i][j] = std::max(prev_job_same_machine_completion, same_job_prev_machine_completion) + p_times[current_job_idx][j];
            }
        }
    }

    double max_lb_for_machine = 0.0;

    for (int j = 0; j < m_machines; ++j) { // Dla każdej maszyny j
        // Czas zakończenia ostatniego uszeregowanego zadania na maszynie j
        double head_time = scheduled_perm.empty() ? 0.0 : C_partial.back()[j];

        // Suma czasów przetwarzania pozostałych zadań na maszynie j
        double sum_remaining_on_machine_j = 0.0;
        for (int job_idx = 0; job_idx < n_total_jobs; ++job_idx) {
            if (!is_scheduled[job_idx]) {
                sum_remaining_on_machine_j += p_times[job_idx][j];
            }
        }

        // Minimalny "ogon" dla nieuszeregowanych zadań na maszynach od j+1 do m
        double min_tail_sum = std::numeric_limits<double>::max();
        bool any_unscheduled = false;
        for (int job_idx = 0; job_idx < n_total_jobs; ++job_idx) {
            if (!is_scheduled[job_idx]) {
                any_unscheduled = true;
                double current_tail_sum = 0.0;
                for (int k = j + 1; k < m_machines; ++k) {
                    current_tail_sum += p_times[job_idx][k];
                }
                min_tail_sum = std::min(min_tail_sum, current_tail_sum);
            }
        }
        if (!any_unscheduled) { // Jeśli wszystkie zadania są uszeregowane
            min_tail_sum = 0.0;
        }
        if (min_tail_sum == std::numeric_limits<double>::max() && any_unscheduled) { // Jeśli j jest ostatnią maszyną, ogon jest 0
            min_tail_sum = 0.0;
        }


        double lb_for_machine_j = head_time + sum_remaining_on_machine_j + min_tail_sum;
        max_lb_for_machine = std::max(max_lb_for_machine, lb_for_machine_j);
    }
    return max_lb_for_machine;
}


// Rekurencyjna funkcja BnB
void bnb_recursive(
        const ProcessingTimes& p_times,
        Permutation& current_perm,
        std::vector<bool>& is_scheduled,
        double& global_upper_bound, // Best Cmax found so far
        Permutation& best_perm_global,
        int n_total_jobs)
{
    // Jeśli current_perm zawiera wszystkie zadania, jest to liść
    if (current_perm.size() == (size_t)n_total_jobs) {
        double cmax_leaf = calculate_cmax(p_times, current_perm);
        if (cmax_leaf < global_upper_bound) {
            global_upper_bound = cmax_leaf;
            best_perm_global = current_perm;
        }
        return;
    }

    // Oblicz dolne ograniczenie dla bieżącego węzła
    double lb = calculate_lower_bound(p_times, current_perm, is_scheduled);

    // Przycinanie (pruning)
    if (lb >= global_upper_bound) {
        return;
    }

    // Rozgałęzienie (Branching)
    for (int job_idx = 0; job_idx < n_total_jobs; ++job_idx) {
        if (!is_scheduled[job_idx]) {
            current_perm.push_back(job_idx);
            is_scheduled[job_idx] = true;

            bnb_recursive(p_times, current_perm, is_scheduled, global_upper_bound, best_perm_global, n_total_jobs);

            // Backtrack
            is_scheduled[job_idx] = false;
            current_perm.pop_back();
        }
    }
}

Permutation branch_and_bound_algorithm(const ProcessingTimes& p_times, double& best_cmax_val) {
    int n_jobs = p_times.size();
    if (n_jobs == 0) {
        best_cmax_val = 0.0;
        return {};
    }



    // Inicjalizacja górnego ograniczenia (UB) za pomocą heurystyki NEH
    Permutation initial_perm_neh = neh_algorithm(p_times, best_cmax_val); // best_cmax_val jest tu UB
    double global_upper_bound = best_cmax_val;
    Permutation best_perm_global = initial_perm_neh;

    Permutation current_perm;
    std::vector<bool> is_scheduled(n_jobs, false);

    bnb_recursive(p_times, current_perm, is_scheduled, global_upper_bound, best_perm_global, n_jobs);

    best_cmax_val = global_upper_bound; // Ostateczny najlepszy Cmax
    return best_perm_global;
}

int main() {


    ProcessingTimes p_times_3x3 = {
            {5, 8, 2},
            {6, 3, 7},
            {4, 9, 5}
    };

    // Dane dla Johnsona (3 zadania, 2 maszyny)
    ProcessingTimes p_times_3x2 = {
            //   M1, M2
            {5, 2}, // Zadanie 1
            {1, 6}, // Zadanie 2
            {9, 7}  // Zadanie 3
    };

    double cmax_val;
    Permutation result_perm;

    std::cout << "--- Dane 3x3 ---\n";
    result_perm = brute_force_algorithm(p_times_3x3, cmax_val);
    print_solution("Przeglad Zupelny (3x3)", result_perm, cmax_val);

    result_perm = neh_algorithm(p_times_3x3, cmax_val);
    print_solution("NEH/FNEH (3x3)", result_perm, cmax_val);

    result_perm = fneh_algorithm(p_times_3x3, cmax_val);
    print_solution("FNEH (3x3)", result_perm, cmax_val);

    result_perm = branch_and_bound_algorithm(p_times_3x3, cmax_val);
    print_solution("Branch and Bound (3x3)", result_perm, cmax_val);

    std::cout << "--- Dane 3x2 dla Johnsona ---\n";
    result_perm = johnson_algorithm(p_times_3x2, cmax_val);
    print_solution("Johnson (3x2)", result_perm, cmax_val);

    // Przykładowe dane dla większej liczby zadań (np. 4 zadania, 3 maszyny)
    // Tylko dla heurystyk i BnB (Brute force będzie wolny)
    ProcessingTimes p_times_4x3 = {
            {21, 53, 80}, // Z1
            {65, 10, 23}, // Z2
            {43, 71, 39}, // Z3
            {29, 30, 62}  // Z4
    };

    std::cout << "--- Dane 4x3 ---\n";
    // result_perm = brute_force_algorithm(p_times_4x3, cmax_val); // Może być wolne
    // print_solution("Przeglad Zupelny (4x3)", result_perm, cmax_val);

    result_perm = neh_algorithm(p_times_4x3, cmax_val);
    print_solution("NEH/FNEH (4x3)", result_perm, cmax_val);

    result_perm = fneh_algorithm(p_times_4x3, cmax_val);
    print_solution("FNEH (4x3)", result_perm, cmax_val);

    result_perm = branch_and_bound_algorithm(p_times_4x3, cmax_val);
    print_solution("Branch and Bound (4x3)", result_perm, cmax_val);

    return 0;
}