#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

// Szablon funkcji do wyznaczania współczynników (działa dla float, double, long double)
template <typename T>
std::vector<T> get_wilkinson_coeffs(int n) {
    std::vector<T> coeffs = {1.0};
    for (int r = 1; r <= n; ++r) {
        std::vector<T> new_coeffs(coeffs.size() + 1, 0.0);
        for (size_t i = 0; i < coeffs.size(); ++i) {
            new_coeffs[i] += coeffs[i] * static_cast<T>(-r);
            new_coeffs[i + 1] += coeffs[i];
        }
        coeffs = new_coeffs;
    }
    return coeffs;
}

// Szablon dla schematu Hornera
template <typename T>
T horner(const std::vector<T>& coeffs, T x) {
    T result = 0.0;
    // Idziemy od końca tablicy (od najwyższej potęgi)
    for (int i = coeffs.size() - 1; i >= 0; --i) {
        result = result * x + coeffs[i];
    }
    return result;
}

// Funkcja testująca dla konkretnego typu
template <typename T>
void run_test(const std::string& type_name) {
    int N = 20;
    std::vector<T> a = get_wilkinson_coeffs<T>(N);
    std::vector<T> test_points = {1.0, 20.0}; // Punkty z treści zadania
    
    std::cout << "--- Typ: " << type_name << " ---" << std::endl;
    for (T x : test_points) {
        T value_horner = horner(a, x);
        std::cout << "Dla x = " << std::setw(2) << x 
                  << " | Wartość (błąd): " << std::scientific << value_horner << std::endl;
    }
    
    std::cout << std::endl;
}

int main() {
    // Ustawienie formatowania wyjścia, żeby łatwo było czytać duże liczby
    std::cout << std::setprecision(6);
    
    // Uruchamiamy test dla trzech typów wymaganych w zadaniu
    run_test<float>("float (32-bit, pojedyncza precyzja)");
    run_test<double>("double (64-bit, podwójna precyzja - to samo co Python)");
    run_test<long double>("long double (80-bit lub 128-bit, rozszerzona precyzja)");
    
    return 0;
}