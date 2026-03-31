# Contributing to LISFLOOD-FP

Thank you for your interest in contributing to LISFLOOD-FP! This document provides guidelines and instructions for contributing to the project.

## Code of Conduct

Please be respectful and considerate of others when contributing to this project. All contributors are expected to adhere to professional conduct.

## How to Contribute

### Reporting Bugs

1. Ensure the bug hasn't already been reported by searching GitHub Issues
2. If you're unable to find an open issue addressing the problem, open a new one with:
   - A clear title and description
   - As much relevant information as possible
   - A code sample or test case demonstrating the issue
   - The expected behavior
   - The actual behavior

### Suggesting Enhancements

1. Open a new GitHub Issue with:
   - A clear title and description
   - The rationale for the enhancement
   - How it would benefit other users
   - Any design or implementation ideas you have

### Pull Requests

1. Fork the repository
2. Create a new branch (`git checkout -b feature/your-feature-name`)
3. Make your changes
4. Run the tests to ensure your changes don't break existing functionality
5. Commit your changes with clear, descriptive commit messages
6. Push your branch to your fork
7. Open a Pull Request against the main branch

## Development Guidelines

### Coding Standards

The LISFLOOD-FP codebase follows these coding conventions:

1. **Naming Conventions**:
   - camelCase for variables and functions
   - PascalCase for classes and types
   - ALL_CAPS for constants and macros

2. **Memory Management**:
   - Always check pointers for nullptr before deletion
   - Use smart pointers for new code where appropriate
   - Follow RAII principles for resource management
   - Ensure proper cleanup in all code paths

3. **Error Handling**:
   - Include clear error messages
   - Handle edge cases gracefully
   - Check function return values

4. **Platform Independence**:
   - Abstract platform-specific code
   - Use conditional compilation with feature detection
   - Avoid hardcoded paths or assumptions

### Documentation

When contributing code, please:

1. Document all public functions, classes, and significant code blocks
2. Update relevant documentation files if you're changing behavior
3. Include examples for new features
4. Keep documentation up-to-date with code changes

### Testing

1. Make sure your code passes all existing tests
2. Add new tests for new functionality
3. Consider edge cases in your tests
4. For performance-critical code, include benchmark tests

## License

By contributing to LISFLOOD-FP, you agree that your contributions will be licensed under the same [GNU GPL License](LICENSE) that covers the project.