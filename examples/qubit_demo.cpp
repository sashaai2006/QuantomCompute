#include "../include/qubit.hpp"
#include <iostream>
#include <vector>
#include <map>
#include <iomanip>

using namespace quantum;

/**
 * @brief Демонстрация основных операций с кубитами
 */
void demonstrate_basic_operations() {
    std::cout << "=== Основные операции с кубитами ===\n\n";
    
    // Создание базисных состояний
    auto q0 = Qubit::zero();
    auto q1 = Qubit::one();
    auto qplus = Qubit::plus();
    auto qminus = Qubit::minus();
    
    std::cout << "1. Базисные состояния:\n";
    std::cout << "   |0⟩ = " << q0 << "\n";
    std::cout << "   |1⟩ = " << q1 << "\n";
    std::cout << "   |+⟩ = " << qplus << "\n";
    std::cout << "   |-⟩ = " << qminus << "\n\n";
    
    // Вероятности измерения
    std::cout << "2. Вероятности измерения:\n";
    std::cout << "   P(|+⟩ → 0) = " << qplus.prob_zero() << "\n";
    std::cout << "   P(|+⟩ → 1) = " << qplus.prob_one() << "\n";
    std::cout << "   P(|-⟩ → 0) = " << qminus.prob_zero() << "\n";
    std::cout << "   P(|-⟩ → 1) = " << qminus.prob_one() << "\n\n";
}

/**
 * @brief Демонстрация квантовых гейтов
 */
void demonstrate_quantum_gates() {
    std::cout << "=== Квантовые гейты ===\n\n";
    
    // Гейты Паули
    std::cout << "1. Гейты Паули:\n";
    auto q = Qubit::zero();
    std::cout << "   Начальное состояние: " << q << "\n";
    
    q.apply_x();
    std::cout << "   После X-гейта: " << q << "\n";
    
    q.apply_y();
    std::cout << "   После Y-гейта: " << q << "\n";
    
    q.apply_z();
    std::cout << "   После Z-гейта: " << q << "\n\n";
    
    // Гейт Адамара
    std::cout << "2. Гейт Адамара:\n";
    q = Qubit::zero();
    std::cout << "   |0⟩ → H → " << q.apply_hadamard() << "\n";
    
    q = Qubit::one();
    std::cout << "   |1⟩ → H → " << q.apply_hadamard() << "\n\n";
    
    // Повороты
    std::cout << "3. Повороты на сфере Блоха:\n";
    q = Qubit::zero();
    std::cout << "   Поворот на π/4 вокруг Y: " << q.rotate_y(M_PI/4) << "\n";
    
    q = Qubit::zero();
    std::cout << "   Поворот на π/2 вокруг X: " << q.rotate_x(M_PI/2) << "\n\n";
}

/**
 * @brief Статистический анализ измерений
 */
void demonstrate_measurement_statistics() {
    std::cout << "=== Статистика измерений ===\n\n";
    
    const int num_measurements = 10000;
    
    // Измерение состояния |+⟩
    std::map<int, int> results;
    for (int i = 0; i < num_measurements; ++i) {
        auto q = Qubit::plus();
        results[q.measure()]++;
    }
    
    std::cout << "1. Измерение |+⟩ (" << num_measurements << " раз):" << std::endl;
    std::cout << "   Частота |0⟩: " << results[0] << " (" 
              << std::fixed << std::setprecision(3) 
              << 100.0 * results[0] / num_measurements << "%)" << std::endl;
    std::cout << "   Частота |1⟩: " << results[1] << " ("
              << 100.0 * results[1] / num_measurements << "%)" << std::endl;
    std::cout << "   Ожидаемые частоты: 50% / 50%\n\n";
    
    // Состояние с коэффициентами 0.8|0⟩ + 0.6|1⟩
    results.clear();
    for (int i = 0; i < num_measurements; ++i) {
        auto q = Qubit(Complex(0.8, 0), Complex(0.6, 0));
        results[q.measure()]++;
    }
    
    std::cout << "2. Измерение 0.8|0⟩ + 0.6|1⟩:" << std::endl;
    std::cout << "   Частота |0⟩: " << results[0] << " ("
              << 100.0 * results[0] / num_measurements << "%)" << std::endl;
    std::cout << "   Частота |1⟩: " << results[1] << " ("
              << 100.0 * results[1] / num_measurements << "%)" << std::endl;
    std::cout << "   Ожидаемые частоты: 64% / 36%\n\n";
}

/**
 * @brief Демонстрация алгоритма Дейча
 */
void demonstrate_deutsch_algorithm() {
    std::cout << "=== Алгоритм Дейча ===\n\n";
    
    std::cout << "Квантовый алгоритм для определения константности функции\n";
    std::cout << "Классически: 2 вычисления, Квантово: 1 вычисление\n\n";
    
    // Тестируем f(x) = 0 (константная)
    std::cout << "1. Тест константной функции f(x) = 0:\n";
    
    auto q = Qubit::zero().apply_hadamard();  // Суперпозиция
    std::cout << "   После H: " << q << "\n";
    
    // Применяем оракул f(x) = 0 (ничего не делаем)
    
    q.apply_hadamard();  // Финальный Адамар
    std::cout << "   После финального H: " << q << "\n";
    std::cout << "   Результат измерения: " << q.measure() 
              << " (константная функция)\n\n";
    
    // Тестируем f(x) = x (сбалансированная)
    std::cout << "2. Тест сбалансированной функции f(x) = x:\n";
    
    q = Qubit::zero().apply_hadamard();
    std::cout << "   После H: " << q << "\n";
    
    // Применяем оракул f(x) = x (эквивалент X-гейту)
    q.apply_x();
    
    q.apply_hadamard();
    std::cout << "   После финального H: " << q << "\n";
    std::cout << "   Результат измерения: " << q.measure() 
              << " (сбалансированная функция)\n\n";
}

/**
 * @brief Основная функция
 */
int main() {
    std::cout << "Квантовые вычисления: Класс Qubit\n";
    std::cout << "====================================\n\n";
    
    demonstrate_basic_operations();
    demonstrate_quantum_gates();
    demonstrate_measurement_statistics();
    demonstrate_deutsch_algorithm();
    
    std::cout << "Программа завершена.\n";
    return 0;
}