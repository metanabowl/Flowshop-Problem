#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <limits>
#include <iomanip>


using ProcessingTimes = std::vector<std::vector<double>>;
using Permutation = std::vector<int>;

// Funkcja do obliczania Cmax (maksymalnego czasu zakończenia) dla danej permutacji zadań
// Cmax to czas zakończenia ostatniego zadania na ostatniej maszynie
// processing_times[job][machine] oznacza czas wykonania zadania na danej maszynie
// permutation[i] = job oznacza kolejnosc wykonania zadań
double calculate_cmax(const ProcessingTimes& p_times, const Permutation& perm) {
    // Sprawdzenie przypadków brzegowych: pusta permutacja lub brak danych o czasach przetwarzania
    if (perm.empty()) {
        return 0.0;
    }
    if (p_times.empty()) {
        return 0.0;
    }

    int n_jobs_in_perm = perm.size();      // Liczba zadań w bieżącej permutacji
    int m_machines = p_times[0].size();     // Liczba maszyn (zakładamy, że wszystkie zadania mają taką samą liczbę operacji)

    if (m_machines == 0) {
        return 0.0;
    }

    // completion_times[i][j] = czas zakończenia i-tego zadania w permutacji na maszynie j.
    // Inicjalizacja macierzy czasów zakończenia zerami.
    std::vector<std::vector<double>> C(n_jobs_in_perm, std::vector<double>(m_machines, 0.0));

    // Wypełnianie macierzy czasów zakończenia zgodnie ze wzorem rekurencyjnym
    for (int i = 0; i < n_jobs_in_perm; ++i) { // Iteracja po zadaniach w danej permutacji
        int current_job_idx = perm[i];        // Pobranie indeksu oryginalnego zadania z permutacji
        for (int j = 0; j < m_machines; ++j) { // Iteracja po maszynach
            // Czas zakończenia poprzedniego zadania na tej samej maszynie
            double prev_job_same_machine_completion = (i == 0) ? 0.0 : C[i - 1][j];
            // Czas zakończenia tego samego zadania na poprzedniej maszynie
            double same_job_prev_machine_completion = (j == 0) ? 0.0 : C[i][j - 1];

            // Obliczenie czasu zakończenia bieżącej operacji
            C[i][j] = std::max(prev_job_same_machine_completion, same_job_prev_machine_completion) + p_times[current_job_idx][j];
        }
    }
    // Cmax to czas zakończenia ostatniego zadania na ostatniej maszynie
    return C[n_jobs_in_perm - 1][m_machines - 1];
}

// Pomocnicza funkcja do wypisywania wyniku
void print_solution(const std::string& algorithm_name, const Permutation& perm, double cmax) {
    std::cout << algorithm_name << ":\n";
    std::cout << "  Permutacja: ";
    for (size_t i = 0; i < perm.size(); ++i) {
        // Wyświetlamy zadania 1-indeksowane dla lepszej czytelności
        std::cout << perm[i] + 1 << (i == perm.size() - 1 ? "" : " -> ");
    }
    // Formatowanie Cmax do dwóch miejsc po przecinku
    std::cout << "\n  Cmax: " << std::fixed << std::setprecision(2) << cmax << "\n\n";
}

//  Przegląd zupełny wszystkich permutacji
Permutation brute_force_algorithm(const ProcessingTimes& p_times, double& best_cmax_val) {
    int n_jobs = p_times.size();
    if (n_jobs == 0) {
        best_cmax_val = 0.0;
        return {};
    }

    // Inicjalizacja bieżącej permutacji (0, 1, ..., n-1)
    Permutation current_perm(n_jobs);
    std::iota(current_perm.begin(), current_perm.end(), 0); // Wypełnia sekwencją 0, 1, ..., n-1

    Permutation best_perm = current_perm;                          // Inicjalizacja najlepszej permutacji
    best_cmax_val = std::numeric_limits<double>::max();    // Ustawienie początkowego Cmax na maksymalną wartość

    // Iteracja przez wszystkie permutacje
    do {
        double cmax = calculate_cmax(p_times, current_perm); // Oblicz Cmax dla bieżącej permutacji
        if (cmax < best_cmax_val) {
            best_cmax_val = cmax;           // Zaktualizuj najlepszy Cmax
            best_perm = current_perm;       // Zaktualizuj najlepszą permutację
        }
    } while (std::next_permutation(current_perm.begin(), current_perm.end())); // Generuj następną permutację

    return best_perm;
}

// Klasyczny algorytm NEH (heurystyka konstrukcyjna)
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

    // Krok 1: Oblicz sumę czasów przetwarzania dla każdego zadania na wszystkich maszynach.
    // job_total_times przechowuje pary {suma czasów, indeks zadania}.
    std::vector<std::pair<double, int>> job_total_times(n_jobs);
    for (int i = 0; i < n_jobs; ++i) {
        double total_time = 0;
        for (int j = 0; j < m_machines; ++j) {
            total_time += p_times[i][j];
        }
        job_total_times[i] = {total_time, i};
    }

    // Krok 2: Posortuj zadania malejąco według sumy czasów przetwarzania.
    std::sort(job_total_times.rbegin(), job_total_times.rend());

    Permutation current_best_perm; // Aktualna najlepsza permutacja budowana iteracyjnie

    // Krok 3 i 4: Iteracyjnie dodawaj zadania do permutacji
    for (int k = 0; k < n_jobs; ++k) {
        int job_to_insert = job_total_times[k].second; // Zadanie do wstawienia

        Permutation best_insertion_perm;            // Najlepsza permutacja po wstawieniu bieżącego zadania
        double min_cmax_for_insertion = std::numeric_limits<double>::max(); // Minimalny Cmax dla bieżącego kroku

        // Wypróbuj wstawienie zadania na każdą możliwą pozycję (od 0 do k)
        for (int pos = 0; pos <= k; ++pos) {
            Permutation temp_perm = current_best_perm; // Kopia obecnej częściowej permutacji
            temp_perm.insert(temp_perm.begin() + pos, job_to_insert); // Wstaw zadanie na daną pozycję

            double cmax = calculate_cmax(p_times, temp_perm); // Oblicz Cmax dla nowej permutacji

            // Jeśli znaleziono lepszą pozycję wstawienia
            if (cmax < min_cmax_for_insertion) {
                min_cmax_for_insertion = cmax;      // Zaktualizuj minimalny Cmax dla bieżącego kroku
                best_insertion_perm = temp_perm;     // Zaktualizuj najlepszą permutację dla bieżącego kroku
            }
        }
        current_best_perm = best_insertion_perm; // Zaktualizuj główną permutację
        best_cmax_val = min_cmax_for_insertion;  // Zaktualizuj najlepszy Cmax
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
    // Algorytm Johnsona wymaga dokładnie dwóch maszyn
    if (p_times[0].size() != 2) {
        std::cerr << "Algorytm Johnsona wymaga 2 maszyn. Podano " << p_times[0].size() << " maszyn.\n";
        best_cmax_val = std::numeric_limits<double>::max(); // Błąd, więc zwracamy dużą wartość Cmax
        return {};
    }

    // Struktura do przechowywania informacji o zadaniu
    struct JobInfo {
        int id;     // Oryginalny indeks zadania
        double p1, p2; // Czasy przetwarzania na maszynie M1 i M2
    };

    std::vector<JobInfo> jobs(n_jobs);
    for (int i = 0; i < n_jobs; ++i) {
        jobs[i] = {i, p_times[i][0], p_times[i][1]};
    }

    std::vector<int> group1, group2; // Grupy zadań dla sortowania Johnsona

    // Dzielenie zadań na dwie grupy
    for (const auto& job : jobs) {
        if (job.p1 < job.p2) {
            group1.push_back(job.id); // Zadanie idzie do grupy 1, jeśli p1 < p2
        } else {
            group2.push_back(job.id); // Zadanie idzie do grupy 2, jeśli p1 >= p2
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
    // Skonkatenuj posortowane grupy, aby uzyskać optymalną permutację
    optimal_perm.insert(optimal_perm.end(), group1.begin(), group1.end());
    optimal_perm.insert(optimal_perm.end(), group2.begin(), group2.end());

    // Oblicz Cmax dla znalezionej permutacji
    best_cmax_val = calculate_cmax(p_times, optimal_perm);
    return optimal_perm;
}

//  FNEH - ulepszony NEH z lokalną optymalizacją po każdym kroku
Permutation fneh_algorithm(const ProcessingTimes& p_times, double& best_cmax_val) {
    int n_jobs = p_times.size();
    int m_machines = p_times[0].size();
    if (n_jobs == 0 || m_machines == 0) {
        best_cmax_val = 0.0;
        return {};
    }

    // Krok 1: Oblicz sumy czasów przetwarzania dla każdego zadania
    std::vector<std::pair<double, int>> job_total_times(n_jobs);
    for (int i = 0; i < n_jobs; ++i) {
        job_total_times[i].first = std::accumulate(p_times[i].begin(), p_times[i].end(), 0.0);
        job_total_times[i].second = i;
    }

    // Krok 2: Posortuj zadania malejąco według sumy czasów
    std::sort(job_total_times.rbegin(), job_total_times.rend());

    Permutation perm; // Bieżąca najlepsza permutacja
    best_cmax_val = 0.0; // Inicjalizacja Cmax

    for (int k = 0; k < n_jobs; ++k) {
        int job_to_add = job_total_times[k].second; // Zadanie do dodania w bieżącej iteracji

        Permutation best_partial_perm;             // Najlepsza częściowa permutacja po wstawieniu
        double best_partial_cmax = std::numeric_limits<double>::max(); // Najlepszy Cmax dla bieżącej iteracji

        // Spróbuj dodać nowe zadanie na każdą możliwą pozycję
        for (int pos = 0; pos <= perm.size(); ++pos) {
            Permutation temp = perm; // Kopia bieżącej permutacji
            temp.insert(temp.begin() + pos, job_to_add); // Wstaw zadanie
            double cmax = calculate_cmax(p_times, temp); // Oblicz Cmax

            if (cmax < best_partial_cmax) {
                best_partial_cmax = cmax;
                best_partial_perm = temp;
            }
        }

        perm = best_partial_perm;       // Zaktualizuj permutację
        best_cmax_val = best_partial_cmax; // Zaktualizuj Cmax

        // Akceleracja FNEH: usuń ostatnio dodane zadanie i przetestuj ponownie jego wstawienie
        // Jest to optymalizacja, która próbuje poprawić znalezione rozwiązanie w danym kroku
        int last_added_job = job_to_add;
        // Usuń ostatnio dodane zadanie z permutacji
        perm.erase(std::find(perm.begin(), perm.end(), last_added_job));

        Permutation fneh_best_perm;
        double fneh_best_cmax = std::numeric_limits<double>::max();

        // Ponownie wstaw ostatnio dodane zadanie na najlepszą pozycję w zmodyfikowanej permutacji
        for (int pos = 0; pos <= perm.size(); ++pos) {
            Permutation temp = perm;
            temp.insert(temp.begin() + pos, last_added_job);
            double cmax = calculate_cmax(p_times, temp);
            if (cmax < fneh_best_cmax) {
                fneh_best_cmax = cmax;
                fneh_best_perm = temp;
            }
        }

        perm = fneh_best_perm;       // Zaktualizuj permutację po akceleracji
        best_cmax_val = fneh_best_cmax; // Zaktualizuj Cmax po akceleracji
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
    // To jest podobne do calculate_cmax, ale dla częściowej permutacji
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

    double max_lb_for_machine = 0.0; // Największe dolne ograniczenie spośród wszystkich maszyn

    // Dla każdej maszyny oblicz dolne ograniczenie
    for (int j = 0; j < m_machines; ++j) {
        // Czas zakończenia ostatniego uszeregowanego zadania na maszynie j (head_time)
        double head_time = scheduled_perm.empty() ? 0.0 : C_partial.back()[j];

        // Suma czasów przetwarzania pozostałych (nieuszeregowanych) zadań na maszynie j
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
                for (int k = j + 1; k < m_machines; ++k) { // Sumowanie czasów na kolejnych maszynach
                    current_tail_sum += p_times[job_idx][k];
                }
                min_tail_sum = std::min(min_tail_sum, current_tail_sum);
            }
        }
        // Obsługa przypadku, gdy wszystkie zadania są już uszeregowane
        if (!any_unscheduled) {
            min_tail_sum = 0.0;
        }
        // Obsługa przypadku, gdy j jest ostatnią maszyną (ogon jest 0)
        if (min_tail_sum == std::numeric_limits<double>::max() && any_unscheduled) {
            min_tail_sum = 0.0;
        }

        // Dolne ograniczenie dla bieżącej maszyny
        double lb_for_machine_j = head_time + sum_remaining_on_machine_j + min_tail_sum;
        // Wybieramy największe LB spośród wszystkich maszyn jako ogólne dolne ograniczenie
        max_lb_for_machine = std::max(max_lb_for_machine, lb_for_machine_j);
    }
    return max_lb_for_machine;
}

// Rekurencyjna funkcja BnB
void bnb_recursive(
        const ProcessingTimes& p_times,
        Permutation& current_perm,
        std::vector<bool>& is_scheduled,
        double& global_upper_bound, // Najlepszy Cmax znaleziony do tej pory (górne ograniczenie)
        Permutation& best_perm_global,
        int n_total_jobs) {
    // Warunek bazowy: Jeśli current_perm zawiera wszystkie zadania, jest to pełna permutacja (liść drzewa)
    if (current_perm.size() == (size_t)n_total_jobs) {
        double cmax_leaf = calculate_cmax(p_times, current_perm); // Oblicz Cmax dla tej permutacji
        if (cmax_leaf < global_upper_bound) {
            global_upper_bound = cmax_leaf;    // Zaktualizuj globalne górne ograniczenie
            best_perm_global = current_perm;    // Zaktualizuj globalnie najlepszą permutację
        }
        return; // Zakończ rekurencję dla tej gałęzi
    }

    // Oblicz dolne ograniczenie (Lower Bound - LB) dla bieżącego węzła (częściowej permutacji)
    double lb = calculate_lower_bound(p_times, current_perm, is_scheduled);

    // Przycinanie (pruning): Jeśli dolne ograniczenie jest większe lub równe
    // aktualnemu górnemu ograniczeniu, to ta gałąź nie może doprowadzić do lepszego rozwiązania.
    if (lb >= global_upper_bound) {
        return; // Nie ma sensu dalej rozwijać tej gałęzi
    }

    // Rozgałęzienie (Branching): Próbuj dodać każde nieuszeregowane zadanie do bieżącej permutacji
    for (int job_idx = 0; job_idx < n_total_jobs; ++job_idx) {
        if (!is_scheduled[job_idx]) { // Jeśli zadanie nie zostało jeszcze uszeregowane
            current_perm.push_back(job_idx);    // Dodaj zadanie do bieżącej permutacji
            is_scheduled[job_idx] = true;       // Oznacz zadanie jako uszeregowane

            // Wywołanie rekurencyjne dla nowego stanu (głębszego węzła)
            bnb_recursive(p_times, current_perm, is_scheduled, global_upper_bound, best_perm_global, n_total_jobs);

            // Backtrack: Cofnij zmiany, aby móc eksplorować inne gałęzie
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

    // Krok 1: Inicjalizacja górnego ograniczenia (UB) za pomocą heurystyki (np. NEH)
    // Dobre początkowe górne ograniczenie może znacznie przyspieszyć działanie BnB.
    Permutation initial_perm_neh = neh_algorithm(p_times, best_cmax_val); // best_cmax_val zostaje tu ustawione jako UB
    double global_upper_bound = best_cmax_val;
    Permutation best_perm_global = initial_perm_neh; // Najlepsza permutacja znaleziona do tej pory (początkowo z NEH)

    Permutation current_perm;                    // Bieżąca częściowa permutacja
    std::vector<bool> is_scheduled(n_jobs, false); // Wektor do śledzenia uszeregowanych zadań

    // Krok 2: Wywołanie rekurencyjnej funkcji BnB
    bnb_recursive(p_times, current_perm, is_scheduled, global_upper_bound, best_perm_global, n_jobs);

    best_cmax_val = global_upper_bound; // Ostateczny, optymalny Cmax
    return best_perm_global;
}

int main() {
    // Przykładowe dane dla 3 zadań i 3 maszyn
    ProcessingTimes p_times_3x3 = {
            {5, 8, 2}, // Czasy dla Zadania 1 na M1, M2, M3
            {6, 3, 7}, // Czasy dla Zadania 2 na M1, M2, M3
            {4, 9, 5}  // Czasy dla Zadania 3 na M1, M2, M3
    };

    // Przykładowe dane dla Johnsona (3 zadania, 2 maszyny)
    ProcessingTimes p_times_3x2 = {
            //   M1, M2
            {5, 2}, // Zadanie 1
            {1, 6}, // Zadanie 2
            {9, 7}  // Zadanie 3
    };

    double cmax_val;         // Zmienna do przechowywania obliczonego Cmax
    Permutation result_perm; // Zmienna do przechowywania znalezionej permutacji

    std::cout << "--- Dane 3x3 ---\n";
    // Testowanie algorytmu Przeglądu Zupełnego
    result_perm = brute_force_algorithm(p_times_3x3, cmax_val);
    print_solution("Przeglad Zupelny (3x3)", result_perm, cmax_val);

    // Testowanie algorytmu NEH
    result_perm = neh_algorithm(p_times_3x3, cmax_val);
    print_solution("NEH (3x3)", result_perm, cmax_val);

    // Testowanie algorytmu FNEH
    result_perm = fneh_algorithm(p_times_3x3, cmax_val);
    print_solution("FNEH (3x3)", result_perm, cmax_val);

    // Testowanie algorytmu Branch and Bound
    result_perm = branch_and_bound_algorithm(p_times_3x3, cmax_val);
    print_solution("Branch and Bound (3x3)", result_perm, cmax_val);

    std::cout << "--- Dane 3x2 dla Johnsona ---\n";
    // Testowanie algorytmu Johnsona (tylko dla 2 maszyn)
    result_perm = johnson_algorithm(p_times_3x2, cmax_val);
    print_solution("Johnson (3x2)", result_perm, cmax_val);

    // Przykładowe dane dla większej liczby zadań (np. 4 zadania, 3 maszyny)
    // Brute force będzie zbyt wolny dla większej liczby zadań
    ProcessingTimes p_times_4x3 = {
            {21, 53, 80}, // Zadanie 1
            {65, 10, 23}, // Zadanie 2
            {43, 71, 39}, // Zadanie 3
            {29, 30, 62}  // Zadanie 4
    };

    std::cout << "--- Dane 4x3 ---\n";
    // brute_force_algorithm jest zakomentowany, ponieważ dla 4 zadań może być już bardzo wolny (4! = 24 permutacji, ale złożoność obliczania Cmax rośnie)
    // result_perm = brute_force_algorithm(p_times_4x3, cmax_val);
    // print_solution("Przeglad Zupelny (4x3)", result_perm, cmax_val);

    // Testowanie NEH dla 4x3
    result_perm = neh_algorithm(p_times_4x3, cmax_val);
    print_solution("NEH (4x3)", result_perm, cmax_val);

    // Testowanie FNEH dla 4x3
    result_perm = fneh_algorithm(p_times_4x3, cmax_val);
    print_solution("FNEH (4x3)", result_perm, cmax_val);

    // Testowanie Branch and Bound dla 4x3
    result_perm = branch_and_bound_algorithm(p_times_4x3, cmax_val);
    print_solution("Branch and Bound (4x3)", result_perm, cmax_val);

    return 0;
}