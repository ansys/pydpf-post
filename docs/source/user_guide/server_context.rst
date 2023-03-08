.. _user_guide_server_context:

==============
Server context
==============

When using a DPF Server released after Anys 2023 R1, a distinction is introduced between two
types of interaction with Ansys licenses:

- When the **Entry** DPF Server context is active, underlying operations will check if a valid Ansys
  license exists, but are not allowed to check-out any license. This means that operations requiring
   to check-out a license are not be available and will raise an error.
- When the **Premium** DPF Server context is active, underlying operations will check if a valid
  Ansys license exists, and are allowed to check-out this license if needed. This means that all DPF
  features are available, but a license may be checked-out.

By default, using PyDPF-Post will start a DPF Server with the **Entry** context active.
To learn more, see the `PyDPF-Core documentation <https://dpf.docs.pyansys.com/dev/user_guide/server_context.html>`_.

Change the default server context
---------------------------------

The default context for the server is **Entry**. You can change the context using
the ``ANSYS_DPF_SERVER_CONTEXT`` environment variable. For more information, see
the `ServerContext class documentation <https://dpf.docs.pyansys.com/dev/api/ansys.dpf.core.server_context.html>`_.
You can also change the server context with this code:

.. code-block::

    from ansys.dpf import post
    post.set_default_server_context(post.AvailableServerContexts.premium)


Release history
---------------

The **Entry** server context is available in server version 6.0
(Ansys 2023 R2) and later.

With a server version earlier than 6.0 (Ansys 2023 R1 and earlier),
**Premium** is the default server context and all features are available,
depending only on their release date.