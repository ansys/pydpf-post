.. _user_guide_server_context:

=========================
Change the server context
=========================

When using a DPF server version later than Ansys 2023 R1, the DPF server context controls Ansys license
interactions.

- When the **Premium** DPF Server context is active, underlying operations check if a valid
  Ansys license exists and are allowed to check out this license if needed. This means that all DPF
  features are available, but a license may be checked-out.
- When the **Entry** DPF Server context is active, underlying operations check if a valid Ansys
  license exists but are not allowed to check out any license. This means that operations that require
  the checkout of a license are not available and raise an error.

By default, using PyDPF-Post starts a DPF server with the **Premium** context active.
For more information, search the `PyDPF-Core documentation <https://dpf.docs.pyansys.com/dev/user_guide/server_context.html>`_
for  **Premium** and **Entry**.

Change the default server context
---------------------------------

The default context for the DPF server is **Premium**. You can change the context using
the ``ANSYS_DPF_SERVER_CONTEXT`` environment variable. For more information, see
`ServerContext <https://dpf.docs.pyansys.com/version/stable/api/ansys.dpf.core.server_context.html>`_ in
the DPF-Core API documentation.

You can also change the DPF server context with this code:

.. code-block::

    from ansys.dpf import post
    post.set_default_server_context(post.AvailableServerContexts.entry)


Release history
---------------

The **Entry** DPF server context is available in DPF Server version 6.0
(Ansys 2023 R2) and later.

With a DPF Server version earlier than 6.0 (Ansys 2023 R1 and earlier),
**Premium** is the default context. All operations are available,
depending only on their release dates.