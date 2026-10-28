# Contribute

Overall guidance on contributing to a PyAnsys library appears in the
[Contributing] topic in the *PyAnsys developer's guide*. Ensure that you
are thoroughly familiar with this guide before attempting to contribute to
PyDPF-Post.

The following contribution information is specific to PyDPF-Post.

[Contributing]: https://dev.docs.pyansys.com/how-to/contributing.html

## Development checks

Run `pre-commit run --all-files` before opening a pull request. Ruff checks
docstrings in place of the archived pydocstyle hook, with the existing rule
exclusions and variadic-parameter behavior preserved. Black and isort continue
to handle formatting and import sorting.

Reusable workflows use the `$/` self-repository syntax. Keep their concurrency
groups distinct from their callers so cancellation does not stop the parent
workflow. Documentation and example runs cancel superseded runs for the same
caller and pull request or ref. Ansys Lab deployments share one concurrency
group because they all update the same destination branch; the latest run wins.
