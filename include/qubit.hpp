#pragma once

#include <complex>
#include <array>
#include <cmath>
#include <random>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace quantum {

using Complex = std::complex<double>;
using StateVector = std::array<Complex, 2>;

/**
 * @brief Класс представляющий квантовый бит (кубит)
 * 
 * Математически кубит представлен как нормированный вектор в ℂ²:
 * |ψ⟩ = α|0⟩ + β|1⟩, где |α|² + |β|² = 1
 * 
 * Информационно: кубит содержит log₂(2^n) = n классических битов информации
 * при измерении, но может находиться в суперпозиции 2^n состояний
 */
class Qubit {
private:
    StateVector state_;  // [α, β] - амплитуды базисных состояний |0⟩ и |1⟩
    
    static thread_local std::mt19937 rng_;  // Генератор случайных чисел для измерений
    static thread_local std::uniform_real_distribution<double> uniform_dist_;

public:
    // Конструкторы
    
    /**
     * @brief Создает кубит в состоянии |0⟩
     */
    Qubit() : state_{Complex(1.0, 0.0), Complex(0.0, 0.0)} {}
    
    /**
     * @brief Создает кубит с заданными амплитудами
     * @param alpha Амплитуда состояния |0⟩
     * @param beta Амплитуда состояния |1⟩
     * @throws std::invalid_argument если состояние не нормировано
     */
    Qubit(const Complex& alpha, const Complex& beta) : state_{alpha, beta} {
        if (!is_normalized()) {
            normalize();
        }
    }
    
    /**
     * @brief Создает кубит из углов на сфере Блоха
     * @param theta Полярный угол θ ∈ [0, π]
     * @param phi Азимутальный угол φ ∈ [0, 2π)
     * 
     * Параметризация: |ψ⟩ = cos(θ/2)|0⟩ + e^(iφ)sin(θ/2)|1⟩
     */
    Qubit(double theta, double phi) {
        state_[0] = std::cos(theta / 2.0);
        state_[1] = std::exp(Complex(0, phi)) * std::sin(theta / 2.0);
    }
    
    // Фабричные методы для стандартных состояний
    
    /**
     * @brief Создает кубит в состоянии |0⟩
     */
    static Qubit zero() {
        return Qubit();
    }
    
    /**
     * @brief Создает кубит в состоянии |1⟩
     */
    static Qubit one() {
        return Qubit(Complex(0.0, 0.0), Complex(1.0, 0.0));
    }
    
    /**
     * @brief Создает кубит в состоянии |+⟩ = (|0⟩ + |1⟩)/√2
     */
    static Qubit plus() {
        const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
        return Qubit(Complex(inv_sqrt2, 0.0), Complex(inv_sqrt2, 0.0));
    }
    
    /**
     * @brief Создает кубит в состоянии |-⟩ = (|0⟩ - |1⟩)/√2
     */
    static Qubit minus() {
        const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
        return Qubit(Complex(inv_sqrt2, 0.0), Complex(-inv_sqrt2, 0.0));
    }
    
    /**
     * @brief Создает кубит в состоянии |i⟩ = (|0⟩ + i|1⟩)/√2
     */
    static Qubit i_state() {
        const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
        return Qubit(Complex(inv_sqrt2, 0.0), Complex(0.0, inv_sqrt2));
    }
    
    // Доступ к состоянию
    
    /**
     * @brief Возвращает амплитуду состояния |0⟩
     */
    const Complex& alpha() const { return state_[0]; }
    
    /**
     * @brief Возвращает амплитуду состояния |1⟩
     */
    const Complex& beta() const { return state_[1]; }
    
    /**
     * @brief Возвращает вектор состояния
     */
    const StateVector& state() const { return state_; }
    
    // Вероятности измерения
    
    /**
     * @brief Возвращает вероятность измерить |0⟩
     */
    double prob_zero() const {
        return std::norm(state_[0]);
    }
    
    /**
     * @brief Возвращает вероятность измерить |1⟩
     */
    double prob_one() const {
        return std::norm(state_[1]);
    }
    
    // Операции измерения
    
    /**
     * @brief Выполняет измерение в вычислительном базисе
     * @return 0 или 1 с соответствующими вероятностями
     * 
     * После измерения состояние коллапсирует к |0⟩ или |1⟩
     */
    int measure() {
        double prob_0 = prob_zero();
        double rand_val = uniform_dist_(rng_);
        
        if (rand_val < prob_0) {
            state_[0] = Complex(1.0, 0.0);
            state_[1] = Complex(0.0, 0.0);
            return 0;
        } else {
            state_[0] = Complex(0.0, 0.0);
            state_[1] = Complex(1.0, 0.0);
            return 1;
        }
    }
    
    /**
     * @brief Выполняет неразрушающее измерение (только возвращает результат)
     * @return 0 или 1 с соответствующими вероятностями
     */
    int peek() const {
        double prob_0 = prob_zero();
        double rand_val = uniform_dist_(rng_);
        return (rand_val < prob_0) ? 0 : 1;
    }
    
    // Унитарные преобразования (квантовые гейты)
    
    /**
     * @brief Применяет гейт Паули-X (NOT)
     * Матрица: [[0, 1], [1, 0]]
     */
    Qubit& apply_x() {
        std::swap(state_[0], state_[1]);
        return *this;
    }
    
    /**
     * @brief Применяет гейт Паули-Y
     * Матрица: [[0, -i], [i, 0]]
     */
    Qubit& apply_y() {
        Complex temp = state_[0];
        state_[0] = Complex(0, -1) * state_[1];
        state_[1] = Complex(0, 1) * temp;
        return *this;
    }
    
    /**
     * @brief Применяет гейт Паули-Z
     * Матрица: [[1, 0], [0, -1]]
     */
    Qubit& apply_z() {
        state_[1] = -state_[1];
        return *this;
    }
    
    /**
     * @brief Применяет гейт Адамара
     * Матрица: (1/√2) * [[1, 1], [1, -1]]
     */
    Qubit& apply_hadamard() {
        const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
        Complex new_alpha = inv_sqrt2 * (state_[0] + state_[1]);
        Complex new_beta = inv_sqrt2 * (state_[0] - state_[1]);
        state_[0] = new_alpha;
        state_[1] = new_beta;
        return *this;
    }
    
    /**
     * @brief Применяет фазовый гейт S
     * Матрица: [[1, 0], [0, i]]
     */
    Qubit& apply_s() {
        state_[1] *= Complex(0, 1);
        return *this;
    }
    
    /**
     * @brief Применяет T-гейт (π/8 гейт)
     * Матрица: [[1, 0], [0, e^(iπ/4)]]
     */
    Qubit& apply_t() {
        state_[1] *= std::exp(Complex(0, M_PI / 4.0));
        return *this;
    }
    
    /**
     * @brief Применяет поворот вокруг оси X
     * @param angle Угол поворота в радианах
     */
    Qubit& rotate_x(double angle) {
        double cos_half = std::cos(angle / 2.0);
        double sin_half = std::sin(angle / 2.0);
        Complex new_alpha = cos_half * state_[0] - Complex(0, sin_half) * state_[1];
        Complex new_beta = cos_half * state_[1] - Complex(0, sin_half) * state_[0];
        state_[0] = new_alpha;
        state_[1] = new_beta;
        return *this;
    }
    
    /**
     * @brief Применяет поворот вокруг оси Y
     * @param angle Угол поворота в радианах
     */
    Qubit& rotate_y(double angle) {
        double cos_half = std::cos(angle / 2.0);
        double sin_half = std::sin(angle / 2.0);
        Complex new_alpha = cos_half * state_[0] - sin_half * state_[1];
        Complex new_beta = cos_half * state_[1] + sin_half * state_[0];
        state_[0] = new_alpha;
        state_[1] = new_beta;
        return *this;
    }
    
    /**
     * @brief Применяет поворот вокруг оси Z
     * @param angle Угол поворота в радианах
     */
    Qubit& rotate_z(double angle) {
        Complex phase_factor = std::exp(Complex(0, angle / 2.0));
        state_[0] *= std::conj(phase_factor);
        state_[1] *= phase_factor;
        return *this;
    }
    
    // Утилитные функции
    
    /**
     * @brief Проверяет нормированность состояния
     */
    bool is_normalized(double tolerance = 1e-10) const {
        double norm_sq = std::norm(state_[0]) + std::norm(state_[1]);
        return std::abs(norm_sq - 1.0) < tolerance;
    }
    
    /**
     * @brief Нормирует состояние
     */
    Qubit& normalize() {
        double norm = std::sqrt(std::norm(state_[0]) + std::norm(state_[1]));
        if (norm > 1e-15) {
            state_[0] /= norm;
            state_[1] /= norm;
        }
        return *this;
    }
    
    /**
     * @brief Вычисляет внутреннее произведение с другим кубитом
     * ⟨ψ|φ⟩ = α₁*α₂ + β₁*β₂
     */
    Complex inner_product(const Qubit& other) const {
        return std::conj(state_[0]) * other.state_[0] + 
               std::conj(state_[1]) * other.state_[1];
    }
    
    /**
     * @brief Вычисляет точность (fidelity) с другим кубитом
     * F(ρ,σ) = |⟨ψ|φ⟩|²
     */
    double fidelity(const Qubit& other) const {
        return std::norm(inner_product(other));
    }
    
    // Операторы
    
    /**
     * @brief Оператор равенства (с учетом глобальной фазы)
     */
    bool operator==(const Qubit& other) const {
        return fidelity(other) > 0.999999; // Учет численных погрешностей
    }
    
    bool operator!=(const Qubit& other) const {
        return !(*this == other);
    }
    
    // Вывод состояния
    
    /**
     * @brief Возвращает строковое представление состояния
     */
    std::string to_string() const {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4);
        
        if (std::abs(state_[0].imag()) < 1e-10 && std::abs(state_[1].imag()) < 1e-10) {
            // Действительные коэффициенты
            oss << state_[0].real() << "|0⟩";
            if (state_[1].real() >= 0) oss << " + ";
            else oss << " ";
            oss << state_[1].real() << "|1⟩";
        } else {
            // Комплексные коэффициенты
            oss << "(" << state_[0].real();
            if (state_[0].imag() >= 0) oss << " + ";
            else oss << " ";
            oss << state_[0].imag() << "i)|0⟩";
            
            oss << " + (" << state_[1].real();
            if (state_[1].imag() >= 0) oss << " + ";
            else oss << " ";
            oss << state_[1].imag() << "i)|1⟩";
        }
        
        return oss.str();
    }
    
    friend std::ostream& operator<<(std::ostream& os, const Qubit& qubit) {
        os << qubit.to_string();
        return os;
    }
};

// Статические члены
thread_local std::mt19937 Qubit::rng_{std::random_device{}()};
thread_local std::uniform_real_distribution<double> Qubit::uniform_dist_{0.0, 1.0};

} // namespace quantum