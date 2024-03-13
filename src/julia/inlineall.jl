import MacroTools as MT

isfundef(fdef) = MT.@capture(MT.longdef1(fdef),
                          function (fcall_ | fcall_) body_ end)
inline_fun(ex) = isfundef(ex) ? :( @inline $ex) : ex

inline_all(ex) = MT.postwalk(inline_fun, ex)

macro inlineall(ex)
    esc(inline_all(ex))
end
