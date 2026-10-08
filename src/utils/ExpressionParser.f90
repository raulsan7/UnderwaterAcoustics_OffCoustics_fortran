! ===============================================================================
! MODULE: ExpressionParser
! ===============================================================================
! Defines the compiler and the evaluator of the mathematical expressions that
! the namelists give as text, e.g. a force as a function of time,
! '1000*sin(2*pi*50*t)', or of the frequency, '100/(1+(f/50)**2)'. Fortran
! cannot evaluate a text at run time, so compile_expression turns it into a
! list of operations in postfix order, checking its syntax, and
! evaluate_expression runs that list over a whole vector of values of the
! variable at once.
!
! Grammar, with the usual precedence (** is right-associative and binds tighter
! than the unary minus, so -t**2 means -(t**2)):
!     sum     := product { (+|-) product }
!     product := unary { (*|/) unary }
!     unary   := (+|-) unary | power
!     power   := primary [ ** unary ]
!     primary := number | name | name '(' sum { ',' sum } ')' | '(' sum ')'
! Numbers are written as 12, 0.5, .5, 1.0e-3 or 2d3. The names are the variable
! (t by default; the caller of compile_expression chooses it, e.g. f), the
! constant pi and the functions of FUNCTION_NAMES. Upper and lower case are the
! same, and blanks are ignored.
! -------------------------------------------------------------------------------
! Author  : Raul Sanz Ramirez
! Contact : raul.sanz.ramirez@upm.es
! Affil.  : ETSIAE, Universidad Politecnica de Madrid
! Project : ERC project off-coustics
! Date    : 10/2026
! ===============================================================================
MODULE ExpressionParser

USE Kinds, ONLY: I32, WP, PI

IMPLICIT NONE

PRIVATE
PUBLIC :: Expression_t, compile_expression, evaluate_expression, function_list


! ===============================================================================
! OPERATIONS AND FUNCTIONS
! The code of function k of FUNCTION_NAMES is FUNCTION_BASE + k
! ===============================================================================
INTEGER(I32), PARAMETER :: MAX_TOKENS    = 256_I32    ! [-] Maximum number of operations of an expression
INTEGER(I32), PARAMETER :: OP_NUMBER     = 1_I32      ! [-] Pushes a number
INTEGER(I32), PARAMETER :: OP_VAR        = 2_I32      ! [-] Pushes the variable (t, f, ...)
INTEGER(I32), PARAMETER :: OP_ADD        = 3_I32      ! [-] a + b
INTEGER(I32), PARAMETER :: OP_SUB        = 4_I32      ! [-] a - b
INTEGER(I32), PARAMETER :: OP_MUL        = 5_I32      ! [-] a * b
INTEGER(I32), PARAMETER :: OP_DIV        = 6_I32      ! [-] a / b
INTEGER(I32), PARAMETER :: OP_POW        = 7_I32      ! [-] a ** b
INTEGER(I32), PARAMETER :: OP_NEG        = 8_I32      ! [-] -a
INTEGER(I32), PARAMETER :: FUNCTION_BASE = 100_I32    ! [-] Code of a function minus its position in FUNCTION_NAMES
INTEGER(I32), PARAMETER :: N_FUNCTIONS   = 22_I32     ! [-] Number of functions

CHARACTER(len=8), PARAMETER :: FUNCTION_NAMES(N_FUNCTIONS) = [CHARACTER(len=8) ::                  &
    'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sinh', 'cosh', 'tanh', 'exp', 'log', 'log10',     &
    'sqrt', 'abs', 'step', 'sgn', 'floor', 'atan2', 'min', 'max', 'mod', 'modulo']    ! [-] Names of the functions

INTEGER(I32), PARAMETER :: FUNCTION_ARGS(N_FUNCTIONS) = [                                           &
    1_I32, 1_I32, 1_I32, 1_I32, 1_I32, 1_I32, 1_I32, 1_I32, 1_I32, 1_I32, 1_I32, 1_I32,            &
    1_I32, 1_I32, 1_I32, 1_I32, 1_I32, 2_I32, 2_I32, 2_I32, 2_I32, 2_I32]            ! [-] Arguments of each function

INTEGER(I32), PARAMETER :: F_SIN     = 1_I32
INTEGER(I32), PARAMETER :: F_COS     = 2_I32
INTEGER(I32), PARAMETER :: F_TAN     = 3_I32
INTEGER(I32), PARAMETER :: F_ASIN    = 4_I32
INTEGER(I32), PARAMETER :: F_ACOS    = 5_I32
INTEGER(I32), PARAMETER :: F_ATAN    = 6_I32
INTEGER(I32), PARAMETER :: F_SINH    = 7_I32
INTEGER(I32), PARAMETER :: F_COSH    = 8_I32
INTEGER(I32), PARAMETER :: F_TANH    = 9_I32
INTEGER(I32), PARAMETER :: F_EXP     = 10_I32
INTEGER(I32), PARAMETER :: F_LOG     = 11_I32
INTEGER(I32), PARAMETER :: F_LOG10   = 12_I32
INTEGER(I32), PARAMETER :: F_SQRT    = 13_I32
INTEGER(I32), PARAMETER :: F_ABS     = 14_I32
INTEGER(I32), PARAMETER :: F_STEP    = 15_I32
INTEGER(I32), PARAMETER :: F_SGN     = 16_I32
INTEGER(I32), PARAMETER :: F_FLOOR   = 17_I32
INTEGER(I32), PARAMETER :: F_ATAN2   = 18_I32
INTEGER(I32), PARAMETER :: F_MIN     = 19_I32
INTEGER(I32), PARAMETER :: F_MAX     = 20_I32
INTEGER(I32), PARAMETER :: F_MOD     = 21_I32
INTEGER(I32), PARAMETER :: F_MODULO  = 22_I32

CHARACTER(len=1), PARAMETER :: END_CHAR = achar(0)    ! [-] Returned by next_char at the end of the text


TYPE :: Expression_t
    ! -------------------------------------------------------------------------
    ! A compiled expression: its operations in postfix order, the value of the
    ! numbers among them, and the depth of the stack that evaluating it needs.
    ! -------------------------------------------------------------------------

    INTEGER(I32) :: n                  = 0_I32       ! [-] Number of operations
    INTEGER(I32) :: depth              = 0_I32       ! [-] Maximum stack depth during the evaluation
    INTEGER(I32) :: code(MAX_TOKENS)   = 0_I32       ! [-] Operation codes
    REAL(WP)     :: value(MAX_TOKENS)  = 0.0_WP      ! [-] Value of the OP_NUMBER operations

END TYPE Expression_t


TYPE :: Parser_t
    ! -------------------------------------------------------------------------
    ! State of compile_expression: the text, the name of its variable, the
    ! position being read, the first error found (empty if none), the depth
    ! of the stack and the expression under construction.
    ! -------------------------------------------------------------------------

    CHARACTER(len=:), ALLOCATABLE :: text         ! [-] Text of the expression
    CHARACTER(len=:), ALLOCATABLE :: variable     ! [-] Name of the variable, in lower case
    CHARACTER(len=:), ALLOCATABLE :: error        ! [-] First error found ('' if none)
    INTEGER(I32)                  :: pos = 1_I32  ! [-] Position being read
    INTEGER(I32)                  :: sp  = 0_I32  ! [-] Current depth of the stack
    TYPE(Expression_t)            :: expr         ! [-] Expression under construction

END TYPE Parser_t


CONTAINS


SUBROUTINE compile_expression(text, expr, error_message, variable)
    ! -------------------------------------------------------------------------
    ! Compiles the text of an expression of one variable, t unless the caller
    ! names another one (e.g. 'f' for a frequency). error_message is empty if
    ! the text is a valid expression, and otherwise describes the first error
    ! and where it is, so that the caller can tell the user which input is
    ! wrong. Any other name than the variable, pi and the functions is an
    ! error, and the message says which variable the expression has.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*)             , INTENT(IN)           :: text             ! [-] Text of the expression
    TYPE(Expression_t)           , INTENT(OUT)          :: expr             ! [-] Compiled expression
    CHARACTER(len=:), ALLOCATABLE, INTENT(OUT)          :: error_message    ! [-] '' if the text is valid
    CHARACTER(len=*)             , INTENT(IN), OPTIONAL :: variable         ! [-] Name of the variable ('t' if absent)

    ! Local variables
    TYPE(Parser_t)   :: p
    CHARACTER(len=1) :: c


    p % text     = trim(adjustl(text))
    p % variable = 't'
    if (present(variable)) p % variable = to_lower(trim(adjustl(variable)))
    p % error    = ''

    if (len(p % text) == 0) then
        CALL set_error(p, 'the expression is empty')
    else
        CALL parse_sum(p)
        if (len(p % error) == 0) then
            c = next_char(p)
            if (c /= END_CHAR) CALL set_error(p, "unexpected '" // c // "'")
        end if
    end if

    if (len(p % error) == 0) expr = p % expr
    error_message = p % error

END SUBROUTINE compile_expression


SUBROUTINE evaluate_expression(expr, x, values)
    ! -------------------------------------------------------------------------
    ! Evaluates a compiled expression for every value of its variable at once:
    ! each operation of the postfix list works on whole vectors. Operations
    ! outside their domain (a logarithm of a negative number, a division by
    ! zero) give NaN or Inf, which the caller must check.
    ! -------------------------------------------------------------------------

    TYPE(Expression_t), INTENT(IN)  :: expr         ! [-] Compiled expression
    REAL(WP)          , INTENT(IN)  :: x(:)         ! [-] Values of the variable (t [s], f [Hz], ...)
    REAL(WP)          , INTENT(OUT) :: values(:)    ! [-] Value of the expression at every x (same size as x)

    ! Local variables
    REAL(WP), ALLOCATABLE :: stack(:,:)
    INTEGER(I32)          :: i, sp


    ALLOCATE(stack(size(x), max(expr % depth, 1_I32)))

    sp = 0_I32
    do i = 1, expr % n
        SELECT CASE (expr % code(i))
        CASE (OP_NUMBER)
            sp = sp + 1_I32
            stack(:, sp) = expr % value(i)
        CASE (OP_VAR)
            sp = sp + 1_I32
            stack(:, sp) = x
        CASE (OP_NEG)
            stack(:, sp) = -stack(:, sp)
        CASE (OP_ADD)
            sp = sp - 1_I32
            stack(:, sp) = stack(:, sp) + stack(:, sp+1)
        CASE (OP_SUB)
            sp = sp - 1_I32
            stack(:, sp) = stack(:, sp) - stack(:, sp+1)
        CASE (OP_MUL)
            sp = sp - 1_I32
            stack(:, sp) = stack(:, sp) * stack(:, sp+1)
        CASE (OP_DIV)
            sp = sp - 1_I32
            stack(:, sp) = stack(:, sp) / stack(:, sp+1)
        CASE (OP_POW)
            sp = sp - 1_I32
            stack(:, sp) = stack(:, sp) ** stack(:, sp+1)
        CASE DEFAULT
            if (FUNCTION_ARGS(expr % code(i) - FUNCTION_BASE) == 2_I32) then
                sp = sp - 1_I32
                CALL apply_function(expr % code(i) - FUNCTION_BASE, stack(:, sp), stack(:, sp+1))
            else
                CALL apply_function(expr % code(i) - FUNCTION_BASE, stack(:, sp))
            end if
        END SELECT
    end do

    values = stack(:, 1)

END SUBROUTINE evaluate_expression


SUBROUTINE apply_function(k, a, b)
    ! -------------------------------------------------------------------------
    ! Applies function k of FUNCTION_NAMES to a, in place. Functions of two
    ! arguments also use b. step(x) is 1 for x >= 0 and 0 otherwise, sgn(x) is
    ! -1, 0 or 1, mod follows Fortran's mod (the sign of a) and modulo follows
    ! Fortran's modulo (the sign of b, handy for periodic signals).
    ! -------------------------------------------------------------------------

    INTEGER(I32), INTENT(IN)           :: k       ! [-] Position of the function in FUNCTION_NAMES
    REAL(WP)    , INTENT(INOUT)        :: a(:)    ! [-] First argument, replaced by the result
    REAL(WP)    , INTENT(IN), OPTIONAL :: b(:)    ! [-] Second argument, for the functions of two


    SELECT CASE (k)
    CASE (F_SIN)
        a = sin(a)
    CASE (F_COS)
        a = cos(a)
    CASE (F_TAN)
        a = tan(a)
    CASE (F_ASIN)
        a = asin(a)
    CASE (F_ACOS)
        a = acos(a)
    CASE (F_ATAN)
        a = atan(a)
    CASE (F_SINH)
        a = sinh(a)
    CASE (F_COSH)
        a = cosh(a)
    CASE (F_TANH)
        a = tanh(a)
    CASE (F_EXP)
        a = exp(a)
    CASE (F_LOG)
        a = log(a)
    CASE (F_LOG10)
        a = log10(a)
    CASE (F_SQRT)
        a = sqrt(a)
    CASE (F_ABS)
        a = abs(a)
    CASE (F_STEP)
        a = merge(1.0_WP, 0.0_WP, a >= 0.0_WP)
    CASE (F_SGN)
        a = merge(1.0_WP, 0.0_WP, a > 0.0_WP) - merge(1.0_WP, 0.0_WP, a < 0.0_WP)
    CASE (F_FLOOR)
        a = a - modulo(a, 1.0_WP)
    CASE (F_ATAN2)
        a = atan2(a, b)
    CASE (F_MIN)
        a = min(a, b)
    CASE (F_MAX)
        a = max(a, b)
    CASE (F_MOD)
        a = mod(a, b)
    CASE (F_MODULO)
        a = modulo(a, b)
    END SELECT

END SUBROUTINE apply_function


FUNCTION function_list() RESULT(list)
    ! -------------------------------------------------------------------------
    ! Names of the functions that an expression can use, separated by blanks.
    ! -------------------------------------------------------------------------

    CHARACTER(len=:), ALLOCATABLE :: list    ! [-] Names of the functions

    ! Local variables
    INTEGER(I32) :: k


    list = ''
    do k = 1, N_FUNCTIONS
        list = list // trim(FUNCTION_NAMES(k)) // ' '
    end do
    list = trim(list)

END FUNCTION function_list


! ===============================================================================
! PARSER (recursive descent, one routine per rule of the grammar)
! ===============================================================================
RECURSIVE SUBROUTINE parse_sum(p)
    ! -------------------------------------------------------------------------
    ! sum := product { (+|-) product }
    ! -------------------------------------------------------------------------

    TYPE(Parser_t), INTENT(INOUT) :: p    ! [-] Parser state

    ! Local variables
    CHARACTER(len=1) :: c


    CALL parse_product(p)

    do while (len(p % error) == 0)
        c = next_char(p)
        if (c /= '+' .and. c /= '-') exit

        p % pos = p % pos + 1_I32
        CALL parse_product(p)
        if (len(p % error) == 0) CALL emit(p, merge(OP_ADD, OP_SUB, c == '+'), 0.0_WP)
    end do

END SUBROUTINE parse_sum


RECURSIVE SUBROUTINE parse_product(p)
    ! -------------------------------------------------------------------------
    ! product := unary { (*|/) unary }
    ! -------------------------------------------------------------------------

    TYPE(Parser_t), INTENT(INOUT) :: p    ! [-] Parser state

    ! Local variables
    CHARACTER(len=1) :: c


    CALL parse_unary(p)

    do while (len(p % error) == 0)
        c = next_char(p)
        if (c /= '*' .and. c /= '/') exit

        p % pos = p % pos + 1_I32
        CALL parse_unary(p)
        if (len(p % error) == 0) CALL emit(p, merge(OP_MUL, OP_DIV, c == '*'), 0.0_WP)
    end do

END SUBROUTINE parse_product


RECURSIVE SUBROUTINE parse_unary(p)
    ! -------------------------------------------------------------------------
    ! unary := (+|-) unary | power
    ! -------------------------------------------------------------------------

    TYPE(Parser_t), INTENT(INOUT) :: p    ! [-] Parser state

    ! Local variables
    CHARACTER(len=1) :: c


    c = next_char(p)

    if (c == '-' .or. c == '+') then
        p % pos = p % pos + 1_I32
        CALL parse_unary(p)
        if (c == '-' .and. len(p % error) == 0) CALL emit(p, OP_NEG, 0.0_WP)
    else
        CALL parse_power(p)
    end if

END SUBROUTINE parse_unary


RECURSIVE SUBROUTINE parse_power(p)
    ! -------------------------------------------------------------------------
    ! power := primary [ ** unary ]. The exponent is a unary, so the operator is
    ! right-associative (a**b**c = a**(b**c)) and a**-b is valid.
    ! -------------------------------------------------------------------------

    TYPE(Parser_t), INTENT(INOUT) :: p    ! [-] Parser state

    ! Local variables
    CHARACTER(len=1) :: c


    CALL parse_primary(p)
    if (len(p % error) /= 0) return

    c = next_char(p)
    if (c == '*' .and. p % pos < len(p % text)) then
        if (p % text(p % pos+1:p % pos+1) == '*') then
            p % pos = p % pos + 2_I32
            CALL parse_unary(p)
            if (len(p % error) == 0) CALL emit(p, OP_POW, 0.0_WP)
        end if
    end if

END SUBROUTINE parse_power


RECURSIVE SUBROUTINE parse_primary(p)
    ! -------------------------------------------------------------------------
    ! primary := number | name | name '(' sum { ',' sum } ')' | '(' sum ')'
    ! -------------------------------------------------------------------------

    TYPE(Parser_t), INTENT(INOUT) :: p    ! [-] Parser state

    ! Local variables
    CHARACTER(len=1) :: c


    c = next_char(p)

    if ((c >= '0' .and. c <= '9') .or. c == '.') then
        CALL parse_number(p)
    else if (is_letter(c)) then
        CALL parse_name(p)
    else if (c == '(') then
        p % pos = p % pos + 1_I32
        CALL parse_sum(p)
        if (len(p % error) == 0) then
            if (next_char(p) == ')') then
                p % pos = p % pos + 1_I32
            else
                CALL set_error(p, "')' expected")
            end if
        end if
    else if (c == END_CHAR) then
        CALL set_error(p, 'the expression ends too soon')
    else
        CALL set_error(p, "unexpected '" // c // "'")
    end if

END SUBROUTINE parse_primary


SUBROUTINE parse_number(p)
    ! -------------------------------------------------------------------------
    ! Reads a number: digits, an optional decimal point and digits, and an
    ! optional exponent (e, E, d or D, an optional sign and digits).
    ! -------------------------------------------------------------------------

    TYPE(Parser_t), INTENT(INOUT) :: p    ! [-] Parser state

    ! Local variables
    INTEGER(I32) :: i, j, ios
    REAL(WP)     :: x


    i = p % pos
    do while (i <= len(p % text))
        if (.not. (p % text(i:i) >= '0' .and. p % text(i:i) <= '9')) exit
        i = i + 1_I32
    end do

    if (i <= len(p % text)) then
        if (p % text(i:i) == '.') then
            i = i + 1_I32
            do while (i <= len(p % text))
                if (.not. (p % text(i:i) >= '0' .and. p % text(i:i) <= '9')) exit
                i = i + 1_I32
            end do
        end if
    end if

    ! The exponent belongs to the number only if digits follow it
    if (i <= len(p % text)) then
        if (index('eEdD', p % text(i:i)) > 0) then
            j = i + 1_I32
            if (j <= len(p % text)) then
                if (p % text(j:j) == '+' .or. p % text(j:j) == '-') j = j + 1_I32
            end if
            if (j <= len(p % text)) then
                if (p % text(j:j) >= '0' .and. p % text(j:j) <= '9') then
                    i = j
                    do while (i <= len(p % text))
                        if (.not. (p % text(i:i) >= '0' .and. p % text(i:i) <= '9')) exit
                        i = i + 1_I32
                    end do
                end if
            end if
        end if
    end if

    read(p % text(p % pos:i-1), *, iostat=ios) x
    if (ios /= 0_I32) then
        CALL set_error(p, "'" // p % text(p % pos:i-1) // "' is not a number")
    else
        p % pos = i
        CALL emit(p, OP_NUMBER, x)
    end if

END SUBROUTINE parse_number


RECURSIVE SUBROUTINE parse_name(p)
    ! -------------------------------------------------------------------------
    ! Reads a name: the variable of the expression, the constant pi, or a
    ! function followed by its arguments between parentheses.
    ! -------------------------------------------------------------------------

    TYPE(Parser_t), INTENT(INOUT) :: p    ! [-] Parser state

    ! Local variables
    CHARACTER(len=:), ALLOCATABLE :: name
    CHARACTER(len=1)              :: c
    INTEGER(I32)                  :: i, k, n_args
    CHARACTER(len=12)             :: n_str


    i = p % pos
    do while (i <= len(p % text))
        if (.not. (is_letter(p % text(i:i)) .or. (p % text(i:i) >= '0' .and. p % text(i:i) <= '9'))) exit
        i = i + 1_I32
    end do
    name    = to_lower(p % text(p % pos:i-1))
    p % pos = i

    if (next_char(p) == '(') then
        k = 0_I32
        do i = 1, N_FUNCTIONS
            if (trim(FUNCTION_NAMES(i)) == name) k = i
        end do
        if (k == 0_I32) then
            CALL set_error(p, "unknown function '" // name // "' (the functions are " // function_list() // ")")
            return
        end if

        p % pos = p % pos + 1_I32
        n_args  = 0_I32
        do
            CALL parse_sum(p)
            if (len(p % error) /= 0) return

            n_args = n_args + 1_I32
            c      = next_char(p)
            if (c == ',') then
                p % pos = p % pos + 1_I32
            else if (c == ')') then
                p % pos = p % pos + 1_I32
                exit
            else
                CALL set_error(p, "',' or ')' expected")
                return
            end if
        end do

        if (n_args /= FUNCTION_ARGS(k)) then
            write(n_str, '(I0)') FUNCTION_ARGS(k)
            CALL set_error(p, "function '" // name // "' takes " // trim(n_str) // " argument(s)")
        else
            CALL emit(p, FUNCTION_BASE + k, 0.0_WP)
        end if
    else if (name == p % variable) then
        CALL emit(p, OP_VAR, 0.0_WP)
    else if (name == 'pi') then
        CALL emit(p, OP_NUMBER, PI)
    else
        CALL set_error(p, "unknown name '" // name // "' (use " // p % variable // ", pi or a function)")
    end if

END SUBROUTINE parse_name


SUBROUTINE emit(p, code, val)
    ! -------------------------------------------------------------------------
    ! Appends an operation to the expression and follows the depth of the
    ! stack: numbers and the variable push one value, binary operations and
    ! functions of two arguments pop one, and the rest leave it as it is.
    ! -------------------------------------------------------------------------

    TYPE(Parser_t), INTENT(INOUT) :: p       ! [-] Parser state
    INTEGER(I32)  , INTENT(IN)    :: code    ! [-] Operation code
    REAL(WP)      , INTENT(IN)    :: val     ! [-] Value, for OP_NUMBER


    if (p % expr % n >= MAX_TOKENS) then
        CALL set_error(p, 'the expression is too long')
        return
    end if

    SELECT CASE (code)
    CASE (OP_NUMBER, OP_VAR)
        p % sp = p % sp + 1_I32
    CASE (OP_ADD:OP_POW)
        p % sp = p % sp - 1_I32
    CASE (FUNCTION_BASE+1:)
        p % sp = p % sp - (FUNCTION_ARGS(code - FUNCTION_BASE) - 1_I32)
    END SELECT

    p % expr % depth            = max(p % expr % depth, p % sp)
    p % expr % n                = p % expr % n + 1_I32
    p % expr % code(p % expr % n)  = code
    p % expr % value(p % expr % n) = val

END SUBROUTINE emit


SUBROUTINE set_error(p, message)
    ! -------------------------------------------------------------------------
    ! Stores the first error found, with its position in the text.
    ! -------------------------------------------------------------------------

    TYPE(Parser_t)  , INTENT(INOUT) :: p          ! [-] Parser state
    CHARACTER(len=*), INTENT(IN)    :: message    ! [-] Description of the error

    ! Local variables
    CHARACTER(len=12) :: pos_str


    if (len(p % error) /= 0) return

    write(pos_str, '(I0)') p % pos
    p % error = message // " at position " // trim(pos_str) // " of '" // p % text // "'"

END SUBROUTINE set_error


FUNCTION next_char(p) RESULT(c)
    ! -------------------------------------------------------------------------
    ! Skips blanks and returns the next character without reading it, or
    ! END_CHAR at the end of the text.
    ! -------------------------------------------------------------------------

    TYPE(Parser_t), INTENT(INOUT) :: p    ! [-] Parser state
    CHARACTER(len=1)              :: c    ! [-] Next character


    do while (p % pos <= len(p % text))
        if (p % text(p % pos:p % pos) /= ' ' .and. p % text(p % pos:p % pos) /= achar(9)) exit
        p % pos = p % pos + 1_I32
    end do

    if (p % pos > len(p % text)) then
        c = END_CHAR
    else
        c = p % text(p % pos:p % pos)
    end if

END FUNCTION next_char


PURE FUNCTION is_letter(c) RESULT(is_it)
    ! -------------------------------------------------------------------------
    ! True for a letter or an underscore, which can start a name.
    ! -------------------------------------------------------------------------

    CHARACTER(len=1), INTENT(IN) :: c        ! [-] Character
    LOGICAL                      :: is_it    ! [-] Whether it is a letter or '_'


    is_it = (c >= 'a' .and. c <= 'z') .or. (c >= 'A' .and. c <= 'Z') .or. c == '_'

END FUNCTION is_letter


PURE FUNCTION to_lower(str) RESULT(res)
    ! -------------------------------------------------------------------------
    ! Lower-case copy of a string.
    ! -------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: str    ! [-] String
    CHARACTER(len=len(str))      :: res    ! [-] The string in lower case

    ! Local variables
    INTEGER(I32) :: i


    res = str
    do i = 1, len(str)
        if (str(i:i) >= 'A' .and. str(i:i) <= 'Z') res(i:i) = achar(iachar(str(i:i)) + 32)
    end do

END FUNCTION to_lower

END MODULE ExpressionParser
