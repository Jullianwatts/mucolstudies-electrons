from pyLCIO import IMPL, IOIMPL, EVENT
#event = EVENT.ParticleID()
impl_pid = IMPL.ParticleIDImpl()
ioimpl_pid = IOIMPL.ParticleIDIOImpl()
#print(type(event))
#print(dir(event))
print("IMPL class type:", type(impl_pid))
print("IOIMPL class type:", type(ioimpl_pid))

print("\nIMPL methods:")
print([m for m in dir(impl_pid) if not m.startswith("__")])

print("\nIOIMPL methods:")
print([m for m in dir(ioimpl_pid) if not m.startswith("__")])

